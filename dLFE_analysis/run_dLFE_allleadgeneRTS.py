#!/usr/bin/env python
# coding: utf-8
import gffutils
import os, glob, pandas as pd, numpy as np, matplotlib.pyplot as plt, logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from codon_randomization import SynonymousCodonPermutingRandomization, NucleotidePermutationRandomization
from rnafold_vienna import RNAfold_direct
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from config import config, init_logging, process_id, ensure_dir_exists

# ---------- Global storage ----------
all_species_delta = {
    'Leading_RTS': []
}

# ---------- Utility Functions ----------
def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def parse_gene_info(gff_file):
    gene_info = {}
    stop_codon_pos = {}
    db = gffutils.create_db(gff_file, ":memory:", merge_strategy="merge")
    for gene in db.features_of_type("gene"):
        chrom = gene.chrom
        start, end = gene.start, gene.end
        strand = gene.strand
        gene_info[gene.id] = (start, end, strand, chrom)
        cds_list = list(db.children(gene, featuretype="CDS", order_by="start"))
        if cds_list:
            cds = cds_list[-1] if strand == "+" else cds_list[0]
            stop_codon_pos[gene.id] = cds.end if strand == "+" else cds.start
    return gene_info, stop_codon_pos

def read_fasta(fasta_file):
    genome_dict = {}
    for rec in SeqIO.parse(fasta_file, "fasta"):
        genome_dict[rec.id] = str(rec.seq)
    return genome_dict

def scan_folder(folder):
    fasta_files = glob.glob(os.path.join(folder, "*.fasta"))
    gff_files = glob.glob(os.path.join(folder, "*.gff3"))
    basename2fasta = {os.path.splitext(os.path.basename(f))[0]: f for f in fasta_files}
    basename2gff = {os.path.splitext(os.path.basename(g))[0]: g for g in gff_files}
    common = sorted(set(basename2fasta) & set(basename2gff))
    return [(bn, basename2fasta[bn], basename2gff[bn]) for bn in common]

# ---------- Identify Leading genes for RTS analysis ----------
def identify_leading_genes(gff_file, fasta_file, prefix, gap_threshold):
    _, leading_regions = extract_leading_genes(gff_file, fasta_file, gap_threshold)
    leading_genes = list(leading_regions.keys())
    logging.info(f"{prefix}: Leading gene identification complete! Found {len(leading_genes)} Leading genes for RTS analysis")
    print(f"{prefix}: {len(leading_genes)} Leading genes identified for RTS analysis")
    return leading_genes

def extract_leading_genes(gene_annotation_file, genome_fasta, gap_threshold):
    processed_gene_pairs = set()
    processed_genes = set()
    leading_stop_codon_regions = {}
    gene_info_dict, stop_codon_pos_dict = parse_gene_info(gene_annotation_file)
    genome_seq = read_fasta(genome_fasta)
    rev_comp_cache = {}
    for chrom, seq in genome_seq.items():
        rev_comp_cache[chrom] = reverse_complement(seq)
    if not gene_info_dict or not genome_seq:
        logging.critical(f"Extract failed: {os.path.basename(gene_annotation_file)}")
        return processed_gene_pairs, leading_stop_codon_regions
    prefix = os.path.splitext(os.path.basename(gene_annotation_file))[0]
    interaction_energy_path = f"{prefix}_leading_RTS_interaction_free_energy.output"
    os.makedirs(os.path.dirname(interaction_energy_path) or '.', exist_ok=True)
    try:
        with open(interaction_energy_path, 'a', encoding='utf-8') as energy_file:
            sorted_genes = sorted(gene_info_dict.items(), key=lambda x: (x[1][3], x[1][0]))
            def _proc(gene, chrom):
                nonlocal leading_stop_codon_regions
                if gene not in processed_genes:
                    processed_genes.add(gene)
                    leading_stop_codon_regions = _process_single_gene_rts(
                        gene, gene_info_dict, genome_seq, rev_comp_cache, stop_codon_pos_dict,
                        energy_file, leading_stop_codon_regions
                    )
            for i in range(len(sorted_genes) - 1):
                gene1, (start1, end1, strand1, chromosome1) = sorted_genes[i]
                gene2, (start2, end2, strand2, chromosome2) = sorted_genes[i + 1]
                if chromosome1 != chromosome2:
                    continue
                if strand1 == '+' and strand2 == '+':
                    gap = start2 - end1
                    if gap > gap_threshold:
                        _proc(gene2, chromosome1)
                elif strand1 == '-' and strand2 == '-':
                    gap = start2 - end1
                    if gap > gap_threshold:
                        _proc(gene1, chromosome1)
                elif strand1 == '-' and strand2 == '+':
                    _proc(gene1, chromosome1)
                    _proc(gene2, chromosome2)
                gene_pair = tuple(sorted([gene1, gene2]))
                processed_gene_pairs.add(gene_pair)
    except Exception as e:
        logging.error(f"Leading‑gene RTS process error: {str(e)}", exc_info=True)
    return processed_gene_pairs, leading_stop_codon_regions

def _process_single_gene_rts(gene, gene_info_dict, genome_seq, rev_comp_cache, stop_codon_dict, energy_file, region_dict):
    if gene not in gene_info_dict or gene not in stop_codon_dict:
        return region_dict
    start, end, strand, chrom = gene_info_dict[gene]
    chr_seq = genome_seq.get(chrom, "")
    if not chr_seq:
        return region_dict
    stop_codon = stop_codon_dict[gene]
    if strand == "+":
        cds3_start = max(0, stop_codon - 90)
        cds3_region = chr_seq[cds3_start:stop_codon]
        utr_end = min(len(chr_seq), stop_codon + 90)
        downstream_utr = chr_seq[stop_codon:utr_end]
    else:
        comp_seq = rev_comp_cache[chrom]
        comp_stop = len(comp_seq) - stop_codon + 1
        cds3_start = max(0, comp_stop - 90)
        cds3_region = comp_seq[cds3_start:comp_stop]
        utr_end = min(len(comp_seq), comp_stop + 90)
        downstream_utr = comp_seq[comp_stop:utr_end]
    if gene not in region_dict:
        region_dict[gene] = (cds3_region, downstream_utr)
        energy_file.write(f"{gene}\t{cds3_region}\t{downstream_utr}\n")
    return region_dict

# ---------- Sequence extraction for RTS (CDS 3' + downstream UTR3) ----------
def extract_gene_sequences(gff_file, fasta_file, gene_list, prefix, out_fna_dir):
    if not gene_list:
        logging.warning(f"{prefix}: Leading gene‑list is empty, skipping seq extract!")
        return None, None
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    db = gffutils.create_db(gff_file, ":memory:", merge_strategy="merge")
    ensure_dir_exists(out_fna_dir)
    cds3_out_path = os.path.join(out_fna_dir, f"{prefix}_Leading_RTS_cds3.fna")
    utr3_out_path = os.path.join(out_fna_dir, f"{prefix}_Leading_RTS_utr3.fna")
    with open(cds3_out_path, "w") as cds3_out, open(utr3_out_path, "w") as utr3_out:
        for gene_id in gene_list:
            try:
                gene = db[gene_id]
                cds_feats = list(db.children(gene, featuretype="CDS", order_by="start"))
                if not cds_feats:
                    continue
                full_cds = "".join([
                    str(genome[cds.chrom].seq[cds.start-1:cds.end]) if cds.strand == "+"
                    else str(Seq(genome[cds.chrom].seq[cds.start-1:cds.end]).reverse_complement())
                    for cds in cds_feats
                ])
                if len(full_cds) == 0:
                    continue
                cds3_seq = full_cds[-90:] if len(full_cds) >= 90 else full_cds.ljust(90, 'N')
                if gene.strand == "+":
                    stop_pos = cds_feats[-1].end
                    utr_start = stop_pos
                    utr_end = min(len(genome[gene.chrom].seq), stop_pos + 90)
                    utr3_seq = str(genome[gene.chrom].seq[utr_start:utr_end])
                else:
                    stop_pos = cds_feats[0].start
                    utr_start = max(0, stop_pos - 90)
                    utr_end = stop_pos
                    utr3_seq = str(genome[gene.chrom].seq[utr_start-1:utr_end-1])
                    utr3_seq = str(Seq(utr3_seq).reverse_complement())
                utr3_seq = utr3_seq.ljust(90, 'N')
                SeqIO.write(SeqRecord(Seq(cds3_seq), id=gene_id, description="CDS3_90bp"), cds3_out, "fasta")
                SeqIO.write(SeqRecord(Seq(utr3_seq), id=gene_id, description="UTR3_90bp"), utr3_out, "fasta")
            except Exception as e:
                logging.warning(f"{prefix}: Extracting Leading RTS {gene_id} failed: {str(e)}")
                continue
    logging.info(f"{prefix}: Leading RTS seq extraction complete: {cds3_out_path} | {utr3_out_path}")
    print(f"{prefix}: Leading RTS (CDS3'+UTR3') sequences extracted")
    return cds3_out_path, utr3_out_path

# ---------- ΔLFE calculation for RTS ----------
def calc_deltaLFE(cds_fn, utr_fn, prefix, global_storage, windowWidth, cdsSpan, utrSpan, randomizationDepth, sequenceSamplingFraction):
    if not cds_fn or not utr_fn or not os.path.exists(cds_fn) or not os.path.exists(utr_fn):
        logging.error(f"{prefix}: RTS files missing, skipping deltaLFE calculation")
        return None
    cdsRand = SynonymousCodonPermutingRandomization(geneticCode=1)
    utrRand = NucleotidePermutationRandomization()
    depth, frac = randomizationDepth, sequenceSamplingFraction
    def get_windows_rts(cds3_seq, utr_seq, span=90, win=30):
        total_windows = cdsSpan + utrSpan
        min_len_needed = total_windows + win - 1
        cds_region = cds3_seq[-span:] if len(cds3_seq) >= span else cds3_seq.ljust(span, 'N')
        utr_region = utr_seq[:span + win - 1].ljust(span + win - 1, 'N')
        full_region = (cds_region + utr_region)[:min_len_needed].ljust(min_len_needed, 'N')
        for i in range(total_windows):
            window = full_region[i:i + win].ljust(win, 'N')
            yield window
    def randomize(seq):
        if len(seq) == 0:
            return ""
        if len(seq) % 3 == 0:
            res = cdsRand.randomize(seq)
            return res if isinstance(res, str) else res[2]
        else:
            res = utrRand.randomize(seq)
            return res if isinstance(res, str) else res[2]
    delta = []
    try:
        cds_records = list(SeqIO.parse(cds_fn, 'fasta'))
        utr_records = list(SeqIO.parse(utr_fn, 'fasta'))
        cds_dict = {rec.id: str(rec.seq) for rec in cds_records}
        utr_dict = {rec.id: str(rec.seq) for rec in utr_records}
        common_genes = set(cds_dict.keys()) & set(utr_dict.keys())
        if not common_genes:
            logging.warning(f"{prefix}: RTS CDS/UTR match failed")
            return None
        for n, gene_id in enumerate(sorted(common_genes)):
            if n % frac != 0:
                continue
            cds3_seq, utr_seq = cds_dict[gene_id], utr_dict[gene_id]
            if len(cds3_seq) < 80 or len(utr_seq) < 80:
                continue
            try:
                LFEs = [RNAfold_direct(win) for win in get_windows_rts(cds3_seq, utr_seq, cdsSpan, windowWidth)]
                randLFEs = []
                for _ in range(depth):
                    rand_cds3 = randomize(cds3_seq)
                    rand_utr = randomize(utr_seq)
                    randLFEs.append([RNAfold_direct(win) for win in get_windows_rts(rand_cds3, rand_utr, cdsSpan, windowWidth)])
                delta.append(np.array(LFEs) - np.nanmean(randLFEs, axis=0))
            except Exception as e:
                logging.warning(f"{prefix}: Skip gene {gene_id}, RNAfold error: {str(e)}")
                continue
        if not delta:
            logging.warning(f"{prefix}: No valid delta‑LFE entries after filtering")
            return None
        delta_arr = np.vstack(delta)
        species_delta_mean = np.nanmean(delta_arr, axis=0)
        global_storage['Leading_RTS'].append(species_delta_mean)
        x = np.arange(-cdsSpan, utrSpan)
        pd.DataFrame({
            'Distance_from_stop_codon (nt)': x,
            'Leading_RTS_DeltaLFE_mean (kcal/mol)': species_delta_mean,
            'Leading_RTS_DeltaLFE_std (kcal/mol)': np.nanstd(delta_arr, axis=0),
            'Species': prefix
        }).to_excel(f'{prefix}_Leading_RTS_deltaLFE.xlsx', index=False)
        logging.info(f"{prefix}: Leading RTS ΔLFE calculation complete")
        print(f"{prefix}: Leading RTS ΔLFE calculation completed")
        return species_delta_mean
    except Exception as e:
        logging.error(f"{prefix}: Leading RTS ΔLFE error: {str(e)}", exc_info=True)
        return None

# ---------- Plot combined heatmap‑line figure for RTS ----------
def plot_heatmap_line_combined(global_storage, cdsSpan, utrSpan, heatmap_cmap, heatmap_nbins, vmin, vmax, output_prefix="all_species"):
    delta_data = global_storage.get('Leading_RTS', [])
    if len(delta_data) == 0:
        logging.error("No Leading RTS data for plotting!")
        return
    merged_delta_arr = np.vstack(delta_data)
    global_mean = np.nanmean(merged_delta_arr, axis=0)
    global_std = np.nanstd(merged_delta_arr, axis=0)
    x = np.arange(-cdsSpan, utrSpan)
    delta_mean_clean = np.nan_to_num(global_mean, nan=0)
    line_color = '#2166ac'
    cmap = mcolors.LinearSegmentedColormap.from_list('custom_RdBu', heatmap_cmap, N=heatmap_nbins)
    fig, (ax_heatmap, ax_line) = plt.subplots(2, 1, figsize=(8, 6),
                                              gridspec_kw={'height_ratios': [1, 5]},
                                              sharex=True)
    heatmap_data = delta_mean_clean.reshape(1, -1)
    im = ax_heatmap.imshow(heatmap_data, aspect='auto', cmap=cmap,
                           extent=[x[0], x[-1], 0, 1],
                           vmin=vmin, vmax=vmax)
    cax = inset_axes(ax_heatmap, width="30%", height="30%",
                     loc='upper right',
                     bbox_to_anchor=(0.6, 1.25, 0.3, 0.3),
                     bbox_transform=ax_heatmap.transAxes,
                     borderpad=0)
    cbar = fig.colorbar(im, cax=cax, orientation='horizontal',
                        ticks=np.linspace(vmin, vmax, 3))
    cbar.set_label('ΔLFE (kcal/mol)', fontsize=9, labelpad=5, y=1.2, rotation=0)
    cbar.ax.tick_params(labelsize=8)
    ax_heatmap.set_yticks([])
    ax_heatmap.set_title('Leading RTS ΔLFE Mean Profile (Heatmap)', fontsize=12, pad=20)
    ax_line.plot(x, global_mean, color=line_color, linewidth=2.5,
                 label=f'Leading RTS Mean (n={len(delta_data)} species)')
    ax_line.fill_between(x, global_mean - global_std, global_mean + global_std,
                         color=line_color, alpha=0.25, label='Leading RTS ±1 SD')
    ax_line.axvline(x=0, color='gray', linestyle='--', linewidth=1.5, alpha=0.7, label='Stop codon')
    ax_line.axhline(y=0, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
    ax_line.set_xlabel('Distance from stop codon (nt)', fontsize=12, labelpad=8)
    ax_line.set_ylabel('ΔLFE (kcal/mol)', fontsize=12, labelpad=8)
    ax_line.set_title('Leading RTS ΔLFE Profile with Variation Range', fontsize=14, pad=10)
    ax_line.set_xlim(-90, 90)
    ax_line.set_xticks(np.arange(-90, 91, 30))
    ax_line.legend(loc='lower left', fontsize=10, framealpha=0.9)
    ax_line.grid(True, alpha=0.2, linestyle='-', linewidth=0.5)
    ax_line.tick_params(axis='both', labelsize=10)
    plt.tight_layout()
    png_path = f'{output_prefix}_Leading_RTS_deltaLFE_heatmap_line.png'
    pdf_path = f'{output_prefix}_Leading_RTS_deltaLFE_heatmap_line.pdf'
    plt.savefig(png_path, dpi=400, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close(fig)
    excel_path = f'{output_prefix}_Leading_RTS_deltaLFE_merged.xlsx'
    pd.DataFrame({
        'Distance_from_stop_codon (nt)': x,
        'Global_Leading_RTS_DeltaLFE_mean (kcal/mol)': global_mean,
        'Global_Leading_RTS_DeltaLFE_std (kcal/mol)': global_std,
        'Valid_Species_Count': len(delta_data)
    }).to_excel(excel_path, index=False)
    logging.info(f"Leading RTS results saved!")
    print(f"\n=== Leading RTS Results Saved ===")
    print(f"Species count: {len(delta_data)}")
    print(f"Figures: {png_path} | {pdf_path}")
    print(f"Data: {excel_path}")

# ---------- Main entry ----------
if __name__ == "__main__":
    log_filename = f"{config.tce.TCE_LOG_BASE_NAME}_LeadingRTS_{process_id()}.log"
    init_logging(log_file_name=log_filename)

    # 全部参数读取config.py
    input_folder = config.tce.INPUT_GENOME_FOLDER
    OUTPUT_FNA_DIR = config.tce.OUT_FNA_LEADING_RTS
    LEADING_LONG_GAP = config.tce.LEADING_LONG_GAP
    windowWidth = config.tce.WINDOW_WIDTH
    cdsSpan = config.tce.CDS_SPAN
    utrSpan = config.tce.UTR_SPAN
    randomizationDepth = config.tce.RANDOMIZATION_DEPTH
    sequenceSamplingFraction = config.tce.SEQUENCE_SAMPLING_FRACTION
    HEATMAP_CMAP_COLORS = config.tce.HEATMAP_CMAP_COLORS
    HEATMAP_N_BINS = config.tce.HEATMAP_N_BINS
    HEATMAP_VMIN = config.tce.HEATMAP_VMIN
    HEATMAP_VMAX = config.tce.HEATMAP_VMAX

    all_species_delta['Leading_RTS'].clear()
    file_pairs = scan_folder(input_folder)
    if not file_pairs:
        logging.critical(f"No matched files in {input_folder}!")
        exit(1)
    print(f"Found {len(file_pairs)} species\n")
    for basename, fasta_file, gff_file in file_pairs:
        print(f"\n===== Processing: {basename} =====")
        leading_genes = identify_leading_genes(gff_file, fasta_file, basename, LEADING_LONG_GAP)
        if leading_genes:
            cds3_fn, utr3_fn = extract_gene_sequences(gff_file, fasta_file, leading_genes, basename, OUTPUT_FNA_DIR)
            calc_deltaLFE(cds3_fn, utr3_fn, basename, all_species_delta,
                          windowWidth, cdsSpan, utrSpan, randomizationDepth, sequenceSamplingFraction)
        else:
            print(f"{basename}: No Leading genes found")
        print(f"===== Finished: {basename} =====")
    print(f"\n===== Plotting Combined Figure =====")
    plot_heatmap_line_combined(all_species_delta, cdsSpan, utrSpan,
                               HEATMAP_CMAP_COLORS, HEATMAP_N_BINS, HEATMAP_VMIN, HEATMAP_VMAX,
                               output_prefix="all_species")
    print(f"\n===== All Leading‑RTS Analysis Completed! =====")
