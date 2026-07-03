#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --account= N/A
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --exclusive
#SBATCH --job-name=liftoff_hyp2sia
#SBATCH --output=logs/liftoff_%j.log

mkdir -p "./logs"

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh

conda activate liftoff-env

TARGET_GENOME="../RNA_seq_align/dsia_1.3.fasta"
REFERENCE_GENOME="..//RNA_seq_align/dhyp_1.0-families.fa"
REFERENCE_GFF="../RNA_seq_align/dhyp.rnd3.kw.sort.gff"
UNMAPPED_FEATURES="./sia_unmapped.txt"
INT_DIR="./hyp2sia_int"
OUTPUT_GFF="../RNA_seq_align/dsia_1.3.gff"

mkdir -p "$INT_DIR"

echo "3, 2, 1"

liftoff \
    -g "$REFERENCE_GFF" \
    -o "$OUTPUT_GFF" \
    -u "$UNMAPPED_FEATURES" \
    -copies \
    -p 8 \
    -dir "$INT_DIR" \
    "$TARGET_GENOME" \
    "$REFERENCE_GENOME" 

conda deactivate 

echo "LIFTOFF!!!!"
