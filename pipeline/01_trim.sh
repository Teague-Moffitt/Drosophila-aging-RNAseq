#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --account= N/A
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --job-name=RNA_trim
#SBATCH --array=0-94
#SBATCH --output=./logs/RNA_trim_%A_%a.log

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh

conda activate fastp_env

INPUT_DIR="./raw_reads"
OUTPUT_DIR="./trimmed_reads"

mkdir -p "$OUTPUT_DIR"

R1_FILES=("$INPUT_DIR"/*_1.fq.gz)
R1="${R1_FILES[$SLURM_ARRAY_TASK_ID]}"
R2="${R1/_1.fq.gz/_2.fq.gz}"
SAMPLE=$(basename "$R1" _1.fq.gz)

fastp \
	-i "$R1" \
	-I "$R2" \
	-o "$OUTPUT_DIR/${SAMPLE}_1_trimmed.fq.gz" \
        -O "$OUTPUT_DIR/${SAMPLE}_2_trimmed.fq.gz" \
        --html "$OUTPUT_DIR/${SAMPLE}_fastp.html" \
        --json "$OUTPUT_DIR/${SAMPLE}_fastp.json" \
        --detect_adapter_for_pe \
        --thread "$SLURM_CPUS_PER_TASK"

conda deactivate

echo "${SAMPLE} trimmed"
