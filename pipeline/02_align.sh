#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --account= N/A
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --exclusive
#SBATCH --job-name=RNA_align
#SBATCH --output=./logs/RNA_align_%A_%a.log
#SBATCH --array=0-94

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh

TRIM_DIR="./trimmed_reads"
ALG_DIR="./aligned_reads"
INDEX="./dhyp_index/dhyp_index"

mkdir -p "$ALG_DIR"

R1_FILES=("$TRIM_DIR"/*_1_trimmed.fq.gz)
        R1="${R1_FILES[$SLURM_ARRAY_TASK_ID]}"
        R2="${R1/_1_trimmed.fq.gz/_2_trimmed.fq.gz}"
        SAMPLE=$(basename "$R1" _1_trimmed.fq.gz)

conda activate hisat2_env

hisat2 \
        -x "$INDEX" \
        -1 "$R1" \
        -2 "$R2" \
        -p 32 \
        --dta \
        2> "$ALG_DIR/${SAMPLE}_hisat2.log" \
| samtools sort -@ 32 -o "$ALG_DIR/${SAMPLE}.bam"

samtools index "$ALG_DIR/${SAMPLE}.bam"

conda deactivate 

echo "${SAMPLE} COMPLETE"
