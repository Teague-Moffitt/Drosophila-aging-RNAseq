#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --account= N/A
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --exclusive
#SBATCH --job-name=read_count
#SBATCH --output=./logs/read_counts_%A_%a.log
#SBATCH --array=0-77

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh

ALG_DIR="./aligned_bams"
CNT_DIR="./read_counts"
ANNOTATION="./D.hyp.gtf"

mkdir -p "$CNT_DIR"

FILES=("$ALG_DIR"/*.bam)
FILE="${FILES[$SLURM_ARRAY_TASK_ID]}"
SAMPLE=$(basename "$FILE" .bam | cut -d'_' -f1)

conda activate subread_env

featureCounts \
        -a "$ANNOTATION" \
        -o "${CNT_DIR}/${SAMPLE}_counts.txt" \
        -p \
	-P \
	-B \
	-T 32 \
	"$FILE"

conda deactivate 

echo "${SAMPLE} COMPLETE"



