#!/bin/bash
# Step 4: Taxonomic classification and lifestyle prediction
# Methods: vConTACT2, Majority Rule, geNomad, CAT, VIBRANT, PhaTYP

set -e
INPUT_FASTA="./03_Viral_Discovery/GVD_LMC_vOTUs.fasta"
OUT_DIR="./04_Viral_Taxonomy"

# 1. vConTACT2 [Ref: 45]
vcontact2 --db 'ProkaryoticViralRefSeq201-Merged' --proteins-fp ${OUT_DIR}/proteins.faa \
    --rel-mode 'Diamond' --pcs-mode MCL --vcs-mode ClusterONE \
    --output-dir ${OUT_DIR}/vcontact2_out

# 2. CAT (Taxonomic voting)
CAT contigs -c ${INPUT_FASTA} -d ./databases/CAT_db -t ./databases/CAT_taxonomy -o ${OUT_DIR}/CAT_out

# 3. Lifestyle Prediction [Ref: 49, 50]
# VIBRANT
VIBRANT_run.py -i ${INPUT_FASTA} -f nucl -t 50 -virome -folder ${OUT_DIR}/vibrant_out

# PhaTYP
phabox2 --task phatyp --contigs ${INPUT_FASTA} --outpth ${OUT_DIR}/phatyp_out

# (Note: Use a Python script to consolidate lifestyle results based on consensus)