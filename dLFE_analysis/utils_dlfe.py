# utils_dlfe.py
import os
import glob
import logging
import numpy as np
import pandas as pd
import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import gffutils
from codon_randomization import SynonymousCodonPermutingRandomization,NucleotidePermutationRandomization
from rnafold_vienna import RNAfold_direct
from config import config, ensure_dir_exists

def reverse_complement(seq:str)->str:
    return str(Seq(seq).reverse_complement())

def parse_gene_info(gff_file):
    gene_info = {}
    start_codon_pos = {}
    db = gffutils.create_db(gff_file,":memory:",merge_strategy="merge")
    for gene in db.features_of_type("gene"):
        chrom=gene.chrom
        start,end=gene.start,gene.end
        strand=gene.strand
        gene_info[gene.id]=(start,end,strand,chrom)
        cds_list=list(db.children(gene,featuretype="CDS",order_by="start"))
        if cds_list:
            cds = cds_list[0] if strand=="+" else cds_list[-1]
            start_codon_pos[gene.id]=cds.start if strand=="+" else cds.end
    return gene_info,start_codon_pos

def read_fasta(fasta_file):
    d={}
    for rec in SeqIO.parse(fasta_file,"fasta"):
        d[rec.id]=str(rec.seq)
    return d

def scan_folder(folder,pattern="*.fasta"):
    fasta_files = glob.glob(os.path.join(folder,"*.fasta"))
    gff_files = glob.glob(os.path.join(folder,"*.gff3"))
    bn2fa={os.path.splitext(os.path.basename(f))[0]:f for f in fasta_files}
    bn2gff={os.path.splitext(os.path.basename(g))[0]:g for g in gff_files}
    common=sorted(set(bn2fa)&set(bn2gff))
    return [(bn,bn2fa[bn],bn2gff[bn]) for bn in common]

def extract_gene_sequences(gff_file,fasta_file,gene_list,prefix,gene_type:str):
    if not gene_list:
        logging.warning(f"{prefix}: {gene_type} gene list empty")
        return None,None
    genome=SeqIO.to_dict(SeqIO.parse(fasta_file,"fasta"))
    db=gffutils.create_db(gff_file,":memory:",merge_strategy="merge")
    if gene_type in ("UTNI","Leading"):
        out_dir=config.tce.OUT_FNA_UTNI_LEADING
    elif gene_type=="Overlap":
        out_dir=config.tce.OUT_FNA_TeRe
    else:
        raise ValueError(f"unknown gene_type {gene_type}")
    ensure_dir_exists(out_dir)
    cds_path=os.path.join(out_dir,f"{prefix}_{gene_type}_cds.fna")
    utr5_path=os.path.join(out_dir,f"{prefix}_{gene_type}_5utr.fna")
    cds_rec=[]
    utr_rec=[]
    for gene_id in gene_list:
        try:
            gene=db[gene_id]
            cds_feats=list(db.children(gene,featuretype="CDS",order_by="start"))
            if not cds_feats: continue
            full_cds="".join([
                str(genome[cds.chrom].seq[cds.start-1:cds.end]) if cds.strand=="+"
                else str(Seq(genome[cds.chrom].seq[cds.start-1:cds.end]).reverse_complement())
                for cds in cds_feats
            ])
            if len(full_cds)==0: continue
            if gene.strand=="+":
                utr_start=max(0,gene.start-1-90)
                utr_end=gene.start-1
                utr5_seq=str(genome[gene.chrom].seq[utr_start:utr_end])
            else:
                utr_start=gene.end
                utr_end=min(len(genome[gene.chrom].seq),gene.end+90)
                utr5_seq=str(genome[gene.chrom].seq[utr_start:utr_end])
                utr5_seq=reverse_complement(utr5_seq)
            utr5_seq=utr5_seq.ljust(90,"N")
            cds_rec.append(SeqRecord(Seq(full_cds),id=gene_id,description=""))
            utr_rec.append(SeqRecord(Seq(utr5_seq),id=gene_id,description=""))
        except Exception as e:
            logging.warning(f"{prefix} extract {gene_id} fail", exc_info=True)
    if cds_rec and utr_rec:
        SeqIO.write(cds_rec,cds_path,"fasta")
        SeqIO.write(utr_rec,utr5_path,"fasta")
        logging.info(f"{prefix} {gene_type} seq written {cds_path} {utr5_path}")
        return cds_path,utr5_path
    else:
        logging.warning(f"{prefix}: no valid seq for {gene_type}")
        return None,None

def calc_deltaLFE(cds_fn,utr_fn,prefix,gene_type:str,global_storage:dict):
    if not cds_fn or not utr_fn or not os.path.exists(cds_fn) or not os.path.exists(utr_fn):
        logging.error(f"{prefix}: {gene_type} missing fasta")
        return None
    cdsRand=SynonymousCodonPermutingRandomization(geneticCode=1)
    utrRand=NucleotidePermutationRandomization()
    win_width=config.tce.WINDOW_WIDTH
    utr_span=config.tce.UTR_SPAN
    cds_span=config.tce.CDS_SPAN
    depth=config.tce.RANDOMIZATION_DEPTH
    frac=config.tce.SEQUENCE_SAMPLING_FRACTION

 #  def get_windows(seq,span,win):
 #      region=(seq[-span:] if len(seq)>=span else seq.ljust(span,"N"))+seq[:span+win-1].ljust(span+win-1,"N")
 #      for i in range(2*span+win-win):
 #          yield region[i:i+win]
            
    def get_windows(seq, span, win):
        region = (seq[-span:] if len(seq) >= span else seq.ljust(span, 'N')) + seq[:span + win - 1].ljust(span + win - 1, 'N')
        for i in range(2 * span + win - 1 - win + 1):
            yield region[i:i + win]

    def randomize(seq):
        if len(seq)==0: return ""
        if len(seq)%3==0:
            r=cdsRand.randomize(seq)
            return r if isinstance(r,str) else r[2]
        else:
            r=utrRand.randomize(seq)
            return r if isinstance(r,str) else r[2]
    delta=[]
    try:
        utr_dict={r.id:str(r.seq) for r in SeqIO.parse(utr_fn,"fasta")}
        cds_dict={r.id:str(r.seq) for r in SeqIO.parse(cds_fn,"fasta")}
        common=set(utr_dict.keys())&set(cds_dict.keys())
        if not common:
            logging.warning(f"{prefix}: {gene_type} no common gene id")
            return None
        for n,gid in enumerate(sorted(common)):
            if n%frac!=0: continue
            us,cs=utr_dict[gid],cds_dict[gid]
            if len(us)<80: continue
            full=us+cs
            try:
                lfes=[RNAfold_direct(w) for w in get_windows(full,utr_span,win_width)]
                rlfes=[]
                for _ in range(depth):
                    ru=randomize(us)
                    rc=randomize(cs)
                    rlfes.append([RNAfold_direct(w) for w in get_windows(ru+rc,utr_span,win_width)])
                delta.append(np.array(lfes)-np.nanmean(rlfes,axis=0))
            except Exception as e:
                logging.warning(f"{prefix} skip gene {gid} RNAfold error {e}")
        if not delta:
            logging.warning(f"{prefix}: {gene_type} no valid delta entries")
            return None
        delta_arr=np.vstack(delta)
        species_mean=np.nanmean(delta_arr,axis=0)
        global_storage[gene_type].append(species_mean)
        x=np.arange(-utr_span,cds_span)
        pd.DataFrame({
            "Distance_from_start_codon (nt)":x,
            f"{gene_type}_DeltaLFE_mean (kcal/mol)":species_mean,
            f"{gene_type}_DeltaLFE_std (kcal/mol)":np.nanstd(delta_arr,axis=0),
            "Species":prefix
        }).to_excel(f"{prefix}_{gene_type}_deltaLFE.xlsx",index=False)
        logging.info(f"{prefix}: {gene_type} delta‑LFE finished")
        return species_mean
    except Exception as e:
        logging.error(f"{prefix}: {gene_type} calculation fatal {e}",exc_info=True)
        return None

def plot_heatmap_line_combined(gene_type:str,global_storage:dict,output_prefix="all_species"):
    delta_data=global_storage.get(gene_type,[])
    if len(delta_data)==0:
        logging.error(f"{gene_type}: no plotting data")
        return
    merged=np.vstack(delta_data)
    g_mean=np.nanmean(merged,axis=0)
    g_std=np.nanstd(merged,axis=0)
    utr_span=config.tce.UTR_SPAN
    cds_span=config.tce.CDS_SPAN
    x=np.arange(-utr_span,cds_span)
    clean_mean=np.nan_to_num(g_mean,nan=0)
    if gene_type=="UTNI":
        line_color="#238b45"
    elif gene_type=="Leading":
        line_color="#2166ac"
    elif gene_type=="Overlap":
        line_color="#8856a7"
    else:
        line_color="#444444"
    cmap=mcolors.LinearSegmentedColormap.from_list("custom",config.tce.HEATMAP_CMAP_COLORS,N=config.tce.HEATMAP_N_BINS)
    fig,(ax_heat,ax_line)=plt.subplots(2,1,figsize=(8,6),gridspec_kw={"height_ratios":[1,5]},sharex=True)
    hm_data=clean_mean.reshape(1,-1)
    im=ax_heat.imshow(hm_data,aspect="auto",cmap=cmap,extent=[x[0],x[-1],0,1],
                      vmin=config.tce.HEATMAP_VMIN,vmax=config.tce.HEATMAP_VMAX)
    cax=inset_axes(ax_heat,width="30%",height="30%",loc="upper right",
                   bbox_to_anchor=(0.6,1.25,0.3,0.3),bbox_transform=ax_heat.transAxes,borderpad=0)
    cbar=fig.colorbar(im,cax=cax,orientation="horizontal",ticks=np.linspace(config.tce.HEATMAP_VMIN,config.tce.HEATMAP_VMAX,3))
    cbar.set_label("ΔLFE (kcal/mol)",fontsize=9,labelpad=5,y=1.2,rotation=0)
    cbar.ax.tick_params(labelsize=8)
    ax_heat.set_yticks([])
    ax_heat.set_title(f"{gene_type} ΔLFE Mean Profile (Heatmap)",fontsize=12,pad=20)

    ax_line.plot(x,g_mean,color=line_color,linewidth=2.5,label=f"{gene_type} Mean (n={len(delta_data)} species)")
    ax_line.fill_between(x,g_mean-g_std,g_mean+g_std,color=line_color,alpha=0.25,label=f"{gene_type} ±1 SD")
    ax_line.axvline(0,color="gray",linestyle="--",lw=1.5,alpha=0.7)
    ax_line.axhline(0,color="gray",linestyle="--",lw=1.5,alpha=0.7)
    ax_line.set_xlabel("Distance from start codon (nt)",fontsize=12)
    ax_line.set_ylabel("ΔLFE (kcal/mol)",fontsize=12)
    ax_line.set_title(f"{gene_type} ΔLFE Profile with Variation Range",fontsize=14)
    ax_line.set_xlim(-90,90)
    ax_line.set_xticks(np.arange(-90,91,30))
    ax_line.legend(loc="lower left",fontsize=10,framealpha=0.9)
    ax_line.grid(True,alpha=0.2)
    plt.tight_layout()
    png=f"{output_prefix}_{gene_type}_deltaLFE_heatmap_line.png"
    pdf=f"{output_prefix}_{gene_type}_deltaLFE_heatmap_line.pdf"
    plt.savefig(png,dpi=400,bbox_inches="tight",facecolor="white")
    plt.savefig(pdf,format="pdf",bbox_inches="tight",facecolor="white")
    plt.close(fig)
    excel=f"{output_prefix}_{gene_type}_deltaLFE_merged.xlsx"
    pd.DataFrame({
        "Distance_from_start_codon (nt)":x,
        f"Global_{gene_type}_DeltaLFE_mean (kcal/mol)":g_mean,
        f"Global_{gene_type}_DeltaLFE_std (kcal/mol)":g_std,
        f"Valid_{gene_type}_Species_Count":len(delta_data)
    }).to_excel(excel,index=False)
    logging.info(f"{gene_type} saved {png}, {pdf}, {excel}")