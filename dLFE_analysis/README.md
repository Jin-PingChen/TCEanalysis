# dLFE_UTR_RTS_analysis  
The code for computational implementation of *"Metagenomic mining of translational coupling elements enables programmable gene expression of polycistronic systems in E. coli"*. This pipeline identifies the 5'-UTR/-CDS regions for all UTNI gene, TeRe gene, and Leading gene. It calculates ΔLFE (delta-LFE, RNA folding free energy offset), supports multispecies batch analysis and generates publication-ready visualization outputs.    

**Code origin note**: This pipeline is adapted from the computational method published in *Chemla, Y., et al. A possible universal role for mRNA secondary structure in bacterial translation revealed using a synthetic operon. Nat. Commun. 11, 4827 (2020). https://doi.org/10.1038/s41467-020-19291-1*. We reuse and refactor core algorithm logic; third-party dependency modules (codon_randomization.py, data_helpers.py, genome_model.py, gff.py, odb4.py, mysql_mafold.py, and nucleic_compress.py）are directly imported without modification.  

# Project Structure   

dLFE_UTR_RTS_analysis/  
├── README.md                     # Documentation and usage guide  
├── config.py                     # Global configuration: paths, thresholds, window params, plotting styles  
├── codon_randomization.py        # Synonymous-codon permutation for CDS, nucleotide permutation for UTR  
├── local_cache.py                # SQLite-based persistent key-value cache  
├── data_helpers.py               # Helper functions for data loading and processing  
├── genome_model.py               # Genome structure models and gene interval handling  
├── gff.py                        # GFF/GTF file parsing and gene annotation extraction  
├── odb4.py                       # Database utilities for OrthoDB integration  
├── mysql_mafold.py               # MySQL/MAFold database connection utilities  
├── nucleic_compress.py           # Nucleic acid sequence compression and encoding  
├── rnafold_vienna.py             # ViennaRNA RNAfold wrapper for LFE computation  
├── core_extract_allgene.py       # Core gene identification: UTNI/Leading/TeRe  
├── utils_dlfe.py                 # Utility functions: GFF/FASTA parsing, ΔLFE computation, plotting  
├── run_dLFE_allgene.py           # Main entry: batch multi-species TCE analysis  
├── run_dLFE_allleadgeneRTS.py    # Leading-RTS analysis (stop-codon side 3'CDS-3'UTR)  
└── test_genome/                  # Input directory: paired *.fasta and *.gff3 files  
   
# Environment Setup  

For the convenience of researchers seeking seamless utilization, we have adhered to all environmental configurations as previously established in Michael Peeri's released RTS analysis toolkit (https://github.com/michaelpeeri/rnafold-rts-public).    
