#!/bin/bash
# Step 5: Host assignment, ADS profiling, and Immune Network
# Methods: CRISPR (CRT), tRNA, Homology, PADLOC

set -e
MAG_DIR="./02_Assembly_Binning/dRep_Output/dereplicated_genomes"
VOTU_FASTA="./03_Viral_Discovery/GVD_LMC_vOTUs.fasta"
OUT_DIR="./05_Host_ADS"

# --- Host Assignment ---
# 1. CRISPR Spacers
java -cp CRT.jar crt -minNR 3 ${MAG_DIR}/*.fa > ${OUT_DIR}/all_spacers.txt
# Blast spacers vs vOTUs (97% ID, 90% cov, max 1 mismatch)
blastn -task blastn-short -query ${OUT_DIR}/spacers.fa -db ${VOTU_FASTA} \
    -perc_identity 97 -qcov_hsp_perc 90 -max_target_seqs 1 -outfmt 6 \
    -out ${OUT_DIR}/crispr_matches.txt

# 2. tRNA Matching
aragorn -t ${VOTU_FASTA} -o ${OUT_DIR}/votu_tRNA.txt
# Match tRNA vs MAGs (100% ID)

# 3. Homology
blastn -query ${VOTU_FASTA} -db ${MAG_DIR}/combined_MAGs.fa \
    -perc_identity 80 -evalue 1e-5 -outfmt 6 -out ${OUT_DIR}/homology_matches.txt
# (Filter: >1kb, coverage <= 50% of host contig)

# --- ADS Profiling ---
# 4. PADLOC [Ref: 54]
padloc --fasta ${MAG_DIR}/combined_MAGs.fa --db ./databases/padloc_db \
    --outdir ${OUT_DIR}/padloc_out
# (Post-processing: Remove "other"/"adaptation", filter <70% breadth)