#!/bin/bash
# Step 1: Read-based profiling of ARGs and microbial communities
# Methods: metaWRAP (QC), Diamond (CARD), arg_ranker (Risk), MetaPhlAn (Taxonomy)

set -e

# --- Config ---
INPUT_DIR="./00_RawData"
OUT_DIR="./01_Read_Profiling"
DB_CARD="./databases/card_v3.0.8.dmnd"
THREADS=20

mkdir -p ${OUT_DIR}/QC ${OUT_DIR}/ARG_Profile ${OUT_DIR}/Taxonomy

# Loop through samples
for R1 in ${INPUT_DIR}/*_1.fastq.gz; do
    SAMPLE=$(basename ${R1} _1.fastq.gz)
    R2=${INPUT_DIR}/${SAMPLE}_2.fastq.gz
    
    echo "Processing ${SAMPLE}..."

    # 1. Quality Control (metaWRAP) [Ref: 21]
    # Note: Assumes metaWRAP environment is active
    metawrap read_qc -1 ${R1} -2 ${R2} -o ${OUT_DIR}/QC/${SAMPLE} -t ${THREADS} --skip-bmtagger
    
    CLEAN_1=${OUT_DIR}/QC/${SAMPLE}/final_pure_reads_1.fastq
    CLEAN_2=${OUT_DIR}/QC/${SAMPLE}/final_pure_reads_2.fastq

    # 2. ARG Profiling (Diamond vs CARD) [Ref: 22]
    # Params: query cover 75, id 90, evalue 1e-5
    diamond blastx -d ${DB_CARD} -q ${CLEAN_1} -o ${OUT_DIR}/ARG_Profile/${SAMPLE}_card.m8 \
        --query-cover 75 --id 90 --evalue 1e-5 -k 1 -f 6 --threads ${THREADS}

    # 3. High-risk ARG Assessment (arg_ranker) [Ref: 23]
    arg_ranker -i ${OUT_DIR}/ARG_Profile/${SAMPLE}_card.m8 -o ${OUT_DIR}/ARG_Profile/${SAMPLE}_risk

    # 4. Microbial Taxonomy (MetaPhlAn 4) [Ref: 24]
    metaphlan ${CLEAN_1},${CLEAN_2} --input_type fastq --nproc ${THREADS} \
        -o ${OUT_DIR}/Taxonomy/${SAMPLE}_profile.txt --bowtie2out ${OUT_DIR}/Taxonomy/${SAMPLE}.bowtie2.bz2
done