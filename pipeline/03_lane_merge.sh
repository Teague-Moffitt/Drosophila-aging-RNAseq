#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --account= N/A
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --exclusive
#SBATCH --job-name=bam_merge
#SBATCH --output=./logs/bam_merge_%j.log

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh

ALG_DIR="./TE_aligned_reads"
MERGED_DIR="./TE_merged_bams"

mkdir -p "$MERGED_DIR"

SAMPLES=(
    DBCC046I1
    DBCC046I2
    DBCC046I3
    DBCC046I4
    DBCC046I5
    DBCC046I6
    DBCC046I7
    DBCC046I8
    DBCC046I9
    DBCC046I10
    DBCC046I11
    DBCC046I13
    DBCC046I14
    DBCC046I15
    DBCC046I16
    DBCC046I17
    DBCC046I18
)

for SAMPLE in "${SAMPLES[@]}"; do

    echo "=== Merging: $SAMPLE ==="

    BAMS=("$ALG_DIR"/${SAMPLE}_*.bam)

    if [[ ${#BAMS[@]} -lt 2 ]]; then
        echo "WARNING: Less than 2 BAMs found for $SAMPLE, skipping."
        continue
    fi

    samtools merge \
        -@ 32 \
        -f \
        "$MERGED_DIR/${SAMPLE}_merged.bam" \
        "${BAMS[@]}"

    samtools index "$MERGED_DIR/${SAMPLE}_merged.bam"

    echo "=== Done: $SAMPLE ==="

done

echo "All merges complete: $(date)"
