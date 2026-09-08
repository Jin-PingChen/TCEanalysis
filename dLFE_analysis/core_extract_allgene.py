#!/usr/bin/env python
# coding: utf-8
import gffutils
import os
import logging
import Bio
from Bio import SeqIO
from Bio.Seq import Seq

from config import config


def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def parse_gene_info(gff_file):
    gene_info = {}
    start_codon_pos = {}
    db = gffutils.create_db(gff_file, ":memory:", merge_strategy="merge")
    for gene in db.features_of_type("gene"):
        chrom = gene.chrom
        start, end = gene.start, gene.end
        strand = gene.strand
        gene_info[gene.id] = (start, end, strand, chrom)
        cds_list = list(db.children(gene, featuretype="CDS", order_by="start"))
        if cds_list:
            cds = cds_list[0] if strand == "+" else cds_list[-1]
            start_codon_pos[gene.id] = cds.start if strand == "+" else cds.end
    return gene_info, start_codon_pos


def read_fasta(fasta_file):
    genome_dict = {}
    for rec in SeqIO.parse(fasta_file, "fasta"):
        genome_dict[rec.id] = str(rec.seq)
    return genome_dict


# ===================== UTNI gene identification =====================
def identify_utni_genes(gff_file: str, fasta_file: str, prefix: str):
    gene_info_dict, _ = parse_gene_info(gff_file)
    _ = read_fasta(fasta_file)

    utni_min = config.tce.UTNI_INTERGENIC_MIN
    utni_max = config.tce.UTNI_INTERGENIC_MAX

    processed_gene_pairs = set()
    processed_genes = set()
    utni_genes = []

    sorted_genes = sorted(gene_info_dict.items(), key=lambda item: (item[1][3], item[1][0]))

    for i in range(len(sorted_genes) - 1):
        gene1, (start1, end1, strand1, chromosome1) = sorted_genes[i]
        gene2, (start2, end2, strand2, chromosome2) = sorted_genes[i + 1]

        if chromosome1 != chromosome2 or strand1 != strand2:
            continue

        distance = start2 - end1
        if utni_min < distance < utni_max:
            gene_pair = tuple(sorted([gene1, gene2]))
            processed_gene_pairs.add(gene_pair)

            if strand1 == "+":
                utni_gene = gene2
            else:
                utni_gene = gene1

            if utni_gene in processed_genes:
                continue

            processed_genes.add(utni_gene)
            utni_genes.append(utni_gene)

    return utni_genes


# ===================== Leading gene identification =====================
def identify_leading_genes(gene_annotation_file: str, genome_fasta: str):
    processed_gene_pairs = set()
    processed_genes = set()
    leading_start_codon_regions = {}

    gene_info_dict, start_codon_pos_dict = parse_gene_info(gene_annotation_file)
    genome_seq = read_fasta(genome_fasta)
    gap_threshold = config.tce.LEADING_LONG_GAP

    if not gene_info_dict or not genome_seq:
        return processed_gene_pairs, leading_start_codon_regions

    prefix = os.path.splitext(os.path.basename(gene_annotation_file))[0]
    interaction_energy_path = f"{prefix}_leading_interaction_free_energy.output"
    os.makedirs(os.path.dirname(interaction_energy_path) or '.', exist_ok=True)

    def _process_single_gene(gene, gene_info, genome_seq, start_codon_dict, energy_file, region_dict):
        if gene not in gene_info or gene not in start_codon_dict:
            return region_dict
        start, end, strand, chrom = gene_info[gene]
        chr_seq = genome_seq.get(chrom, "")
        if not chr_seq:
            return region_dict
        start_codon = start_codon_dict[gene]
        if strand == "+":
            region_start = max(0, start_codon - 31)
            region_end = min(len(chr_seq), start_codon + 89)
            rna_region = chr_seq[region_start:start_codon - 1]
            downstream_region = chr_seq[start_codon - 1: region_end]
        else:
            comp_seq = reverse_complement(chr_seq)
            comp_start = len(comp_seq) - start_codon
            region_start = max(0, comp_start - 30)
            region_end = min(len(comp_seq), comp_start + 90)
            rna_region = comp_seq[region_start:comp_start]
            downstream_region = comp_seq[comp_start:region_end]
        if gene not in region_dict:
            region_dict[gene] = (rna_region, downstream_region)
            energy_file.write(f"{gene}\t{rna_region}\t{downstream_region}\n")
        return region_dict

    def _proc(gene, chrom):
        nonlocal leading_start_codon_regions
        if gene not in processed_genes:
            processed_genes.add(gene)
            leading_start_codon_regions = _process_single_gene(
                gene, gene_info_dict, genome_seq, start_codon_pos_dict,
                energy_file, leading_start_codon_regions
            )

    try:
        with open(interaction_energy_path, 'a', encoding='utf-8') as energy_file:
            sorted_genes = sorted(gene_info_dict.items(), key=lambda x: (x[1][3], x[1][0]))
            for i in range(len(sorted_genes) - 1):
                gene1, (start1, end1, strand1, chromosome1) = sorted_genes[i]
                gene2, (start2, end2, strand2, chromosome2) = sorted_genes[i + 1]
                if chromosome1 != chromosome2:
                    continue
                chr_seq = genome_seq.get(chromosome1, "")
                if not chr_seq:
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
        logging.warning(f"identify_leading_genes error: {e}", exc_info=True)
    return processed_gene_pairs, leading_start_codon_regions


# ===================== TeRe Overlap gene identification =====================
def identify_overlap_genes(bed_file: str, gff_file: str, prefix: str):
    import re
    strand_dict = {}
    gene_start_dict = {}
    gene_end_dict = {}
    db = gffutils.create_db(gff_file, ":memory:", merge_strategy="merge")
    for gene in db.features_of_type("gene"):
        gene_id = gene.id
        strand_dict[gene_id] = gene.strand
        gene_start_dict[gene_id] = gene.start
        gene_end_dict[gene_id] = gene.end

    processed_gene_pairs = set()
    processed_genes = set()
    overlap_genes = set()
    gene_id_pattern = re.compile(r'ID=([^;]+)')
    gene_name_pattern = re.compile(r'gene=([^;]+)')

    try:
        with open(bed_file, 'r', encoding="utf-8") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                cols = line.split('\t')
                if len(cols) < 8:
                    continue
                gene1_match = gene_id_pattern.search(cols[3]) or gene_name_pattern.search(cols[3])
                gene2_match = gene_id_pattern.search(cols[7]) or gene_name_pattern.search(cols[7])
                gene1 = gene1_match.group(1).strip('"') if gene1_match else None
                gene2 = gene2_match.group(1).strip('"') if gene2_match else None
                if not gene1 or not gene2 or gene1 == gene2:
                    continue
                s1 = strand_dict.get(gene1)
                s2 = strand_dict.get(gene2)
                if not ((s1 == '+' and s2 == '+') or (s1 == '-' and s2 == '-')):
                    continue
                if gene1 not in gene_start_dict or gene2 not in gene_start_dict:
                    continue
                g1s = gene_start_dict[gene1]
                g1e = gene_end_dict[gene1]
                g2s = gene_start_dict[gene2]
                g2e = gene_end_dict[gene2]

                overlap_start = max(g1s, g2s)
                overlap_end = min(g1e, g2e)
                overlap_len = overlap_end - overlap_start + 1
                if overlap_len <= 0:
                    continue

                pair_key = tuple(sorted([gene1, gene2]))
                if pair_key in processed_gene_pairs:
                    continue
                processed_gene_pairs.add(pair_key)

                if s1 == "+":
                    tere_candidate = gene2
                else:
                    tere_candidate = gene1
                if tere_candidate not in processed_genes:
                    processed_genes.add(tere_candidate)
                    overlap_genes.add(tere_candidate)
    except Exception as e:
        logging.warning(f"identify_overlap_genes error: {e}", exc_info=True)

    stats = {
        'overlap_genes': len(overlap_genes),
        'overlap_pairs': len(processed_gene_pairs),
        'total_genes': len(strand_dict)
    }
    return overlap_genes, stats