# Drosophila-aging-RNAseq
## Overview
This repository contains the bioinformatics pipeline and analysis scripts 
for an exploratory RNA-seq study investigating the genetic basis of sexually 
dimorphic aging and melanization in *Drosophila hypocausta* and *D. siamana*. 
This work was conducted in the Wei Lab at the University of British Columbia.

## Biological Background
*D. hypocausta* and *D. siamana* display distinct sexual dimorphism in both 
lifespan and melanization patterning. This study uses whole-transcriptome 
RNA sequencing spanning the lifespan of the flies to identify the genetic basis driving 
sex-specific aging trajectories and melanization phenotypes.

## Experimental Design
- **Species:** *D. hypocausta* and *D. siamana*
- **Tissues:** Carcass and gonad (testis or ovary)
- **Timepoints:** 3, 5, 10, 20, 40, and 60 days post-eclosion
- **Conditions:** Both sexes sequenced at each timepoint
- Note: Due to library preparation issues, not all timepoints are represented 
for both species and tissues.

## Key Findings
This study identified a candidate gene involved in the melanization phenotype and 
demonstrated systematic differences in gene expression between sexes, 
potentially contributing to overall sexual dimorphism in aging trajectories.

## Pipeline Overview
Raw Reads (FASTQ) ->

01_trim.sh # Quality and adaptor trimming with fastp

02_align.sh # Alignment to reference genome using HISAT2

03_lane_merge.sh # Merge BAM files across multiple sequencing lanes

04_count_reads.sh # Read counting with featureCounts

DEseq_clean.R # Differential expression analysis with DESeq2

### Genome Annotation
`05_liftoff.sh` lifts the *D. hypocausta* genome annotation over to the 
*D. siamana* reference (one of the two haplotypes) using Liftoff, enabling cross-species 
comparative analysis.

## Analysis Scripts

| Script | Description |
|--------|-------------|
| `DEseq_clean.R` | Full DESeq2 differential expression pipeline |
| `Avg_Trajectory_Plot.R` | Average gene expression trajectory line plots across aging timepoints |
| `shared_single_traj.R` | Individual gene trajectory plots for two species |
| `Chromosome_Enrichment.R` | Chromosome-level enrichment analysis of DE genes |
| `species_plot.R` | Cross-species expression comparisons |

## Dependencies

**Pipeline tools**
- fastp
- HISAT2
- featureCounts
- Liftoff

**R packages**
- DESeq2
- ggplot2
- dplyr

## Data Availability
Raw sequencing data is not publicly available. Please contact the Wei Lab 
at UBC for data access inquiries.

## Author
Teague Moffitt
Faculty of Forestry, University of British Columbia
Wei Lab — PI: Dr. Kevin Wei

