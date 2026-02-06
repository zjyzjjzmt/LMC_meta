#!/bin/bash
# Step 2: Assembly, Annotation, and Genome-resolved analysis (MAGs)
# Methods: MEGAHIT, Prodigal, Diamond, metaWRAP, dRep, GTDB-Tk, MetaThermo

set -e

INPUT_DIR="./01_Read_Profiling/QC"
OUT_DIR="./02_Assembly_Binning"
THREADS=50

# Databases
DB_CARD="./databases/card.dmnd"
DB_ICE="./databases/ICEs.dmnd"

mkdir -p ${OUT_DIR}

for SAMPLE_DIR in ${INPUT_DIR}/*; do
    SAMPLE=$(basename ${SAMPLE_DIR})
    R1=${SAMPLE_DIR}/final_pure_reads_1.fastq
    R2=${SAMPLE_DIR}/final_pure_reads_2.fastq

    # 1. Assembly (MEGAHIT) [Ref: 25]
    megahit -1 ${R1} -2 ${R2} --min-contig-len 1000 -m 0.9 -t ${THREADS} -o ${OUT_DIR}/Assembly/${SAMPLE}

    # 2. Contig Annotation (Prodigal + Diamond) [Ref: 26-28]
    CONTIGS=${OUT_DIR}/Assembly/${SAMPLE}/final.contigs.fa
    prodigal -i ${CONTIGS} -p meta -a ${OUT_DIR}/Annotation/${SAMPLE}.faa -d ${OUT_DIR}/Annotation/${SAMPLE}.fna
    
    # ARG & MGE identification (Identity >= 80%, Coverage >= 70%)
    diamond blastp -q ${OUT_DIR}/Annotation/${SAMPLE}.faa -d ${DB_CARD} \
        --id 80 --query-cover 70 --evalue 1e-10 -o ${OUT_DIR}/Annotation/${SAMPLE}_card.txt
    
    diamond blastp -q ${OUT_DIR}/Annotation/${SAMPLE}.faa -d ${DB_ICE} \
        --id 80 --query-cover 70 --evalue 1e-10 -o ${OUT_DIR}/Annotation/${SAMPLE}_ice.txt
    
    # 3. Binning (metaWRAP) [Ref: 21]
    metawrap binning -o ${OUT_DIR}/Binning/${SAMPLE} -t ${THREADS} -a ${CONTIGS} \
        --metabat2 --maxbin2 --concoct ${R1} ${R2}

    # 4. Bin Refinement
    metawrap bin_refinement -o ${OUT_DIR}/Refinement/${SAMPLE} -t ${THREADS} \
        -A ${OUT_DIR}/Binning/${SAMPLE}/metabat2_bins/ \
        -B ${OUT_DIR}/Binning/${SAMPLE}/maxbin2_bins/ \
        -C ${OUT_DIR}/Binning/${SAMPLE}/concoct_bins/ \
        -c 70 -x 10
done

# 5. Dereplication (dRep) [Ref: 33]
# Combine all refined bins into one folder first
dRep dereplicate ${OUT_DIR}/dRep_Output -g ${OUT_DIR}/All_Bins/*.fa \
    -p ${THREADS} -comp 70 -con 10

# 6. Taxonomy (GTDB-Tk) [Ref: 34]
gtdbtk classify_wf --genome_dir ${OUT_DIR}/dRep_Output/dereplicated_genomes \
    --out_dir ${OUT_DIR}/GTDBTk_Output -x fa --cpus ${THREADS}

# 7. Optimal Temperature (MetaThermo) [Ref: 38]
# (Run custom MetaThermo script here)