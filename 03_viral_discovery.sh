#!/bin/bash
# Step 3: Identification and quality control of the viral catalogue (GVD_LMC)
# Methods: VirSorter2, VirFinder, geNomad, CheckV, ClusterGenomes

set -e

# Config
CONTIGS_ALL="./02_Assembly_Binning/All_Contigs_5kb.fasta" # Merged contigs >5kb
OUT_DIR="./03_Viral_Discovery"
THREADS=50

mkdir -p ${OUT_DIR}

# 1. Run Identification Tools (Parallel)
echo "Running VirSorter2..."
virsorter run -i ${CONTIGS_ALL} -w ${OUT_DIR}/vs2_out --include-groups dsDNAphage,ssDNA -j ${THREADS} all

echo "Running geNomad..."
genomad end-to-end --cleanup ${CONTIGS_ALL} ${OUT_DIR}/genomad_out ./databases/genomad_db

echo "Running VirFinder..."
Rscript scripts/run_virfinder.R ${CONTIGS_ALL} ${OUT_DIR}/virfinder_out.tsv

# 2. Filter and Merge Candidates (Custom Python Script)
# Logic: VS2 rules OR VF rules OR geNomad rules
python scripts/filter_viral_contigs.py \
    --vs2 ${OUT_DIR}/vs2_out/final-viral-score.tsv \
    --vf ${OUT_DIR}/virfinder_out.tsv \
    --genomad ${OUT_DIR}/genomad_out/*_summary.tsv \
    --fasta ${CONTIGS_ALL} \
    --out ${OUT_DIR}/merged_viral_candidates.fasta

# 3. Deduplication (CD-HIT)
cd-hit-est -i ${OUT_DIR}/merged_viral_candidates.fasta -o ${OUT_DIR}/viral_nr.fasta \
    -c 1.0 -n 10 -M 16000 -T ${THREADS}

# 4. Quality Control (CheckV) [Ref: 43]
checkv end_to_end ${OUT_DIR}/viral_nr.fasta ${OUT_DIR}/checkv_out -t ${THREADS} \
    -d ./databases/checkv-db-v1.4

# Filter CheckV results: Keep if (viral_gene > 0) OR (host_gene <= 1 & len >= 10kb)
# Extract clean sequences (removing host regions)
python scripts/filter_checkv.py --input ${OUT_DIR}/checkv_out --fasta ${OUT_DIR}/viral_nr.fasta \
    --out ${OUT_DIR}/high_quality_viral.fasta

# 5. Clustering (MIUViG) [Ref: 44]
# 95% ANI, 85% AF
python scripts/ClusterGenomes.py -f ${OUT_DIR}/high_quality_viral.fasta \
    -o ${OUT_DIR}/GVD_LMC_vOTUs.fasta -i 95 -c 85