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

#subset for non DC Mullers
echo -e "MullerA,MullerA\nMullerB,MullerB\nMullerE,MullerE\nMullerF,MullerF" > chroms_nonDC.txt

grep -v "MullerDC" ../dhyp.rnd3.kw.sort.gff > dhyp_nonDC.gff

TARGET_GENOME="../dsia_1.3.fasta"
REFERENCE_GENOME="../dhyp_1.0.fa"
REFERENCE_GFF="dhyp_nonDC.gff"
UNMAPPED_FEATURES="unmapped_nonDC_sia.txt"
INT_DIR="int_nonDC_sia"
OUTPUT_GFF="sia_nonDC.gff"
CHROMS="chroms_nonDC.txt"

mkdir -p "$INT_DIR"

echo "3, 2, 1"

liftoff \
    -g "$REFERENCE_GFF" \
    -o "$OUTPUT_GFF" \
    -u "$UNMAPPED_FEATURES" \
    -copies \
    -p 8 \
    -chroms "$CHROMS" \
    -dir "$INT_DIR" \
    "$TARGET_GENOME" \
    "$REFERENCE_GENOME" 

#subset for DC1/2 haplotypes A and B

grep "MullerDC_1" ../dhyp.rnd3.kw.sort.gff > dhyp_DC1.gff
grep "MullerDC_2" ../dhyp.rnd3.kw.sort.gff > dhyp_DC2.gff

echo "MullerDC_2,MullerDC1.hapA" > chroms_DC1_hapA.txt
echo "MullerDC_2,MullerDC1.hapB" > chroms_DC1_hapB.txt
echo "MullerDC_1,MullerDC2.hapA" > chroms_DC2_hapA.txt
echo "MullerDC_1,MullerDC2.hapB" > chroms_DC2_hapB.txt

#DC1 hapA
TARGET_GENOME="../dsia_1.3.fasta"
REFERENCE_GENOME="../dhyp_1.0.fa"
REFERENCE_GFF="dhyp_DC2.gff"
UNMAPPED_FEATURES="unmapped_DC1_hapA_sia.txt"
INT_DIR="int_DC1_hapA_sia"
OUTPUT_GFF="sia_DC1_hapA.gff"
CHROMS="chroms_DC1_hapA.txt"

mkdir -p "$INT_DIR"

echo "3, 2, 1"

liftoff \
    -g "$REFERENCE_GFF" \
    -o "$OUTPUT_GFF" \
    -u "$UNMAPPED_FEATURES" \
    -p 8 \
    -chroms "$CHROMS" \
    -dir "$INT_DIR" \
    "$TARGET_GENOME" \
    "$REFERENCE_GENOME" 

#DC1 hapB
TARGET_GENOME="../dsia_1.3.fasta"
REFERENCE_GENOME="../dhyp_1.0.fa"
REFERENCE_GFF="dhyp_DC2.gff"
UNMAPPED_FEATURES="unmapped_DC1_hapB_sia.txt"
INT_DIR="int_DC1_hapB_sia"
OUTPUT_GFF="sia_DC1_hapB.gff"
CHROMS="chroms_DC1_hapB.txt"

mkdir -p "$INT_DIR"

echo "3, 2, 1"

liftoff \
    -g "$REFERENCE_GFF" \
    -o "$OUTPUT_GFF" \
    -u "$UNMAPPED_FEATURES" \
    -p 8 \
    -chroms "$CHROMS" \
    -dir "$INT_DIR" \
    "$TARGET_GENOME" \
    "$REFERENCE_GENOME" 

#DC2 hapA
TARGET_GENOME="../dsia_1.3.fasta"
REFERENCE_GENOME="../dhyp_1.0.fa"
REFERENCE_GFF="dhyp_DC1.gff"
UNMAPPED_FEATURES="unmapped_DC2_hapA_sia.txt"
INT_DIR="int_DC2_hapA_sia"
OUTPUT_GFF="sia_DC2_hapA.gff"
CHROMS="chroms_DC2_hapA.txt"

mkdir -p "$INT_DIR"

echo "3, 2, 1"

liftoff \
    -g "$REFERENCE_GFF" \
    -o "$OUTPUT_GFF" \
    -u "$UNMAPPED_FEATURES" \
    -p 8 \
    -chroms "$CHROMS" \
    -dir "$INT_DIR" \
    "$TARGET_GENOME" \
    "$REFERENCE_GENOME" 

#DC2 hap B
TARGET_GENOME="../dsia_1.3.fasta"
REFERENCE_GENOME="../dhyp_1.0.fa"
REFERENCE_GFF="dhyp_DC1.gff"
UNMAPPED_FEATURES="unmapped_DC2_hapB_sia.txt"
INT_DIR="int_DC2_hapB_sia"
OUTPUT_GFF="sia_DC2_hapB.gff"
CHROMS="chroms_DC2_hapB.txt"

mkdir -p "$INT_DIR"

echo "3, 2, 1"

liftoff \
    -g "$REFERENCE_GFF" \
    -o "$OUTPUT_GFF" \
    -u "$UNMAPPED_FEATURES" \
    -p 8 \
    -chroms "$CHROMS" \
    -dir "$INT_DIR" \
    "$TARGET_GENOME" \
    "$REFERENCE_GENOME" 

cat sia_DC1_hapA.gff sia_DC1_hapB.gff sia_DC2_hapA.gff sia_DC2_hapB.gff sia_nonDC.gff > ../dsia_1.3.gff

conda deactivate 

echo "LIFTOFF!!!!"
