#!/usr/bin/env python3
# run_all_tce.py
import logging
import os
from config import config, init_logging, process_id
from utils_dlfe import scan_folder, extract_gene_sequences, calc_deltaLFE, plot_heatmap_line_combined
from core_extract_allgene import identify_utni_genes, identify_leading_genes, identify_overlap_genes

# ===================== Module toggle switches =====================
RUN_UTNI = True
RUN_LEADING = True
RUN_TERE = True
# =================================================================

def scan_tere_folder(folder: str):
    import glob
    bed_files = glob.glob(os.path.join(folder, "*.gene.overlap.bed"))
    pairs = []
    for bed in sorted(bed_files):
        basename = os.path.basename(bed).replace(".gene.overlap.bed", "")
        fasta_path = os.path.join(folder, f"{basename}.fasta")
        gff_path = os.path.join(folder, f"{basename}.gff3")
        if os.path.exists(fasta_path) and os.path.exists(gff_path):
            pairs.append((basename, bed, fasta_path, gff_path))
    return pairs


if __name__ == "__main__":
    log_filename = f"{config.tce.TCE_LOG_BASE_NAME}_{process_id()}.log"
    init_logging(log_file_name=log_filename)

    global_storage = {
        "UTNI": [],
        "Leading": [],
        "Overlap": []
    }
    input_folder = config.tce.INPUT_GENOME_FOLDER
    fa_gff_pairs = scan_folder(input_folder)

    # ---------------------- UTNI analysis ----------------------
    if RUN_UTNI:
        logging.info("===== Starting UTNI gene analysis =====")
        if not fa_gff_pairs:
            logging.warning("No fasta/gff3 pairs for UTNI, skip UTNI module")
        else:
            for prefix, fasta_file, gff_file in fa_gff_pairs:
                logging.info(f"[UTNI] processing species: {prefix}")
                utni_gene_list = identify_utni_genes(gff_file, fasta_file, prefix)
                if utni_gene_list:
                    cds_fn, utr_fn = extract_gene_sequences(gff_file, fasta_file, utni_gene_list, prefix, gene_type="UTNI")
                    calc_deltaLFE(cds_fn, utr_fn, prefix, "UTNI", global_storage)
            plot_heatmap_line_combined("UTNI", global_storage, output_prefix="all_species")
        logging.info("===== UTNI analysis finished =====\n")

    # ---------------------- Leading analysis ----------------------
    if RUN_LEADING:
        logging.info("===== Starting Leading gene analysis =====")
        if not fa_gff_pairs:
            logging.warning("No fasta/gff3 pairs for Leading, skip Leading module")
        else:
            for prefix, fasta_file, gff_file in fa_gff_pairs:
                logging.info(f"[Leading] processing species: {prefix}")
                _, leading_regions = identify_leading_genes(gff_file, fasta_file)
                leading_gene_list = list(leading_regions.keys())
                if leading_gene_list:
                    cds_fn, utr_fn = extract_gene_sequences(gff_file, fasta_file, leading_gene_list, prefix, gene_type="Leading")
                    calc_deltaLFE(cds_fn, utr_fn, prefix, "Leading", global_storage)
            plot_heatmap_line_combined("Leading", global_storage, output_prefix="all_species")
        logging.info("===== Leading analysis finished =====\n")

    # ---------------------- TeRe Overlap analysis ----------------------
    if RUN_TERE:
        logging.info("===== Starting TeRe(Overlap) gene analysis =====")
        tere_triplets = scan_tere_folder(input_folder)
        if not tere_triplets:
            logging.warning("No *.gene.overlap.bed found for TeRe, skip TeRe module")
        else:
            for prefix, bed_file, fasta_file, gff_file in tere_triplets:
                logging.info(f"[TeRe] processing species: {prefix}")
                tere_gene_set, _stats = identify_overlap_genes(bed_file, gff_file, prefix)
                tere_gene_list = list(tere_gene_set)
                if tere_gene_list:
                    cds_fn, utr_fn = extract_gene_sequences(gff_file, fasta_file, tere_gene_list, prefix, gene_type="Overlap")
                    calc_deltaLFE(cds_fn, utr_fn, prefix, "Overlap", global_storage)
            plot_heatmap_line_combined("Overlap", global_storage, output_prefix="all_species")
        logging.info("===== TeRe analysis finished =====\n")

    logging.info("===== All selected TCE?dLFE jobs completed =====")
    print("\nDone! Check log?file, fna output dir and png/pdf/excel outputs.")