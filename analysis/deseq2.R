library(DESeq2)
library(here)
#ALL DATA ----

##make metadata ----
metadata <- read.csv(
  here("Data", "RNAseq samples - Copy of Sheet3.csv"), 
  header = TRUE, 
  row.names = 1,
  colClasses = "factor"
)

metadata$Tissue <- trimws(metadata$Tissue)
metadata$Age..days. <- trimws(metadata$Age..days.)
metadata$Sex <- trimws(metadata$Sex)

metadata$Age..days. <- as.numeric(as.character(
  metadata$Age..days.))

metadata$Sex <- as.factor(metadata$Sex)

metadata$Tissue <- as.factor(metadata$Tissue)

##create list of files ----
count_dir <- here("Data", "raw_read_counts")

files <- list.files(
  path = count_dir,
  pattern = "\\.txt$",
  full.names = TRUE
)


##match samples ----
sample_ids <- sub(
  "_counts\\.txt$",
  "",
  basename(files)
)

metadata_order <- rownames(metadata)

files <- files[
  match(metadata_order, sample_ids)
]


##create count matrix ----
count_list <- lapply(files, function(f) {
  
  df <- read.table(
    f, 
    header = TRUE,
    comment.char = "#",
    stringsAsFactors = FALSE
  )
  
  counts <- df[, c(1, ncol(df))]
  
  sample_name <- sub(
  "_counts\\.txt$",
    "",
    basename(f)
  )
  
  colnames(counts) [2] <- sample_name
  
  counts
})


##merge matrix ----
count_matrix <- Reduce(function(x, y) {
  merge(x, y, by = "Geneid")
}, count_list)

rownames(count_matrix) <- count_matrix$Geneid

count_matrix <- count_matrix[, -1]

count_matrix <- as.matrix(count_matrix)


##make dds object ----
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix, 
  colData = metadata,
  design = ~ Species + Sex + Tissue + Age..days. 
)

dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE)




#ALL KW SAMPLES ----

##create KW metadata ----
KW_samples <- rownames(metadata[grepl("^KW",rownames(metadata)), ])

metadata_KW <- metadata[KW_samples, ]


##subset matrix ----
count_matrix_KW <- count_matrix[, KW_samples]


## make DDS object ----

dds_KW <- DESeqDataSetFromMatrix(
  countData = count_matrix_KW, 
  colData = metadata_KW,
  design = ~ Species + Sex + Tissue + Age..days. 
)

dds_KW <- DESeq(dds_KW)
vsd_KW <- vst(dds_KW, blind = TRUE)


## monotonic ----
inc_kw_age <- subset(results(dds_KW, name = "Age..days."), padj < 0.05 & log2FoldChange > 0)

dec_kw_age <- subset(results(dds_KW, name = "Age..days."), padj < 0.05 & log2FoldChange < 0)

#D. HYPOCASTA CARCASS KW0 ----

##create hypo carcass KW metadata ----
metadata_hypo_carcass <- subset(
  metadata, 
  Tissue == "Carcass" & Species == "D. hypocausta"
)

hypo_carcass_KW_samples <- rownames(
  metadata_hypo_carcass[
    grepl(
      "^KW",
      rownames(metadata_hypo_carcass)
    ), 
  ]
)

metadata_hypo_carcass_KW <- metadata_hypo_carcass[hypo_carcass_KW_samples, ]

metadata_hypo_carcass_KW$Age..days. <- as.numeric(as.character(metadata_hypo_carcass_KW$Age..days.))



##subset matrix ----
count_matrix_hypo_carcass_KW <- count_matrix[, hypo_carcass_KW_samples]



##make DDS object hypocarcass ----
dds_hypo_carcass_KW <- DESeqDataSetFromMatrix(
  countData = count_matrix_hypo_carcass_KW, 
  colData = metadata_hypo_carcass_KW,
  design = ~ Sex * Age..days.
)


dds_hypo_carcass_KW <- DESeq(dds_hypo_carcass_KW)

vsd_hypo_carcass_KW <- vst(dds_hypo_carcass_KW, blind = TRUE)


##Monotonic Genes ----
inc_kw_hypo_car_age <- subset(results(dds_hypo_carcass_KW, name = "Age..days."), padj < 0.05 & log2FoldChange > 0)

dec_kw_hypo_car_age <- subset(results(dds_hypo_carcass_KW, name = "Age..days."), padj < 0.05 & log2FoldChange < 0)

top100_kw_inc <- head(inc_kw_hypo_car_age[order((inc_kw_hypo_car_age$padj)),], n = 100)

top100_kw_dec <- head(dec_kw_hypo_car_age[order((dec_kw_hypo_car_age$padj)),], n = 100)

###view top 100 ----
# top_100_sex_dec <- head(decreasing_sex_ordered, n = 100)
# 
# cat(rownames(top_100_sex_dec), sep = "\n")


#D. SIAMANA CARCASS KW0 ----

##create sia carcass KW metadata ----
metadata_sia_carcass <- subset(
  metadata, 
  Tissue == "Carcass" & Species == "D. siamana"
)

sia_carcass_samples <- rownames(metadata_sia_carcass)

metadata_sia_carcass$Age..days. <- as.numeric(as.character(metadata_sia_carcass$Age..days.))



##subset matrix ----
count_matrix_sia_carcass <- count_matrix[, sia_carcass_samples]



##make DDS object siacarcass ----
dds_sia_carcass <- DESeqDataSetFromMatrix(
  countData = count_matrix_sia_carcass, 
  colData = metadata_sia_carcass,
  design = ~ Sex * Age..days.
)


dds_sia_carcass <- DESeq(dds_sia_carcass)

vsd_sia_carcass <- vst(dds_sia_carcass, blind = TRUE)


##Monotonic Genes ----
inc_sia_car_age <- subset(results(dds_sia_carcass, name = "Age..days."), padj < 0.05 & log2FoldChange > 0)

dec_sia_car_age <- subset(results(dds_sia_carcass, name = "Age..days."), padj < 0.05 & log2FoldChange < 0)

top100_kw_inc <- head(inc_sia_car_age[order((inc_sia_car_age$padj)),], n = 100)

top100_kw_dec <- head(dec_sia_car_age[order((dec_sia_car_age$padj)),], n = 100)

###view top 100 ----
# top_100_sex_dec <- head(decreasing_sex_ordered, n = 100)
# 
# cat(rownames(top_100_sex_dec), sep = "\n")


#REPEATS ----

##create list of files ----
count_dir_TE <- here("Data", "TE_raw_reads")

files <- list.files(
  path = count_dir_TE,
  pattern = "\\.txt$",
  full.names = TRUE
)


##match samples ----
sample_ids <- sub(
  "_counts\\.txt$",
  "",
  basename(files)
)

metadata_order <- rownames(metadata)

files <- files[
  match(metadata_order, sample_ids)
]


##create count matrix ----
count_list <- lapply(files, function(f) {
  
  df <- read.table(
    f, 
    header = TRUE,
    skip = 1,
    comment.char = "",
    sep = "\t",
    check.names = FALSE,
    quote = "",
    row.names = NULL
  )
  
  counts <- df[, c(1, ncol(df))]
  
  sample_name <- sub(
    "_counts\\.txt$",
    "",
    basename(f)
  )
  
  colnames(counts) [2] <- sample_name
  
  counts
})


##merge matrix ----
count_matrix_TE <- Reduce(function(x, y) {
  merge(x, y, by = "Geneid")
}, count_list)

rownames(count_matrix_TE) <- count_matrix_TE$Geneid

count_matrix_TE <- count_matrix_TE[, -1]

count_matrix_TE <- as.matrix(count_matrix_TE)


##subset matrix ----
###for hypocausta ----
metadata_hypo <- subset(metadata, Species == "D. hypocausta")

hypo_samples <- rownames(
  metadata[metadata$Species == "D. hypocausta", ]
)

count_matrix_TE_hypo <- count_matrix_TE[, hypo_samples]

### for hypo carcass ----
hypo_carcass_samples <- rownames(
  metadata_hypo[metadata_hypo$Tissue == "Carcass", ]
)

count_matrix_TE_hypo_carcass <- count_matrix_TE[, hypo_carcass_samples]

### for KW samples ----
hypo_carcass_KW_samples <- rownames(
  metadata_hypo_carcass[grepl("^KW", rownames(metadata_hypo_carcass)), ]
)

count_matrix_TE_hypo_carcass_KW <- count_matrix_TE[, hypo_carcass_KW_samples]

te_counts_filtered <- count_matrix_TE_hypo_carcass_KW[!grepl("rRNA", rownames(count_matrix_TE_hypo_carcass_KW)), ]


##make DEseq object ----
### for all ----
library(DESeq2)

TE_dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_TE, 
  colData = metadata,
  design = ~ Species + Sex + Tissue + Age..days. 
)

TE_dds <- DESeq(TE_dds)
TE_vsd <- varianceStabilizingTransformation(TE_dds, blind = TRUE)

### for hypo ----
TE_dds_hypo <- DESeqDataSetFromMatrix(
  countData = count_matrix_TE_hypo, 
  colData = metadata_hypo,
  design = ~ Sex + Tissue + Age..days. 
)

TE_dds_hypo <- DESeq(TE_dds_hypo)
TE_vsd_hypo <- varianceStabilizingTransformation(TE_dds_hypo, blind = TRUE)

### for hypo carcass ----
TE_dds_hypo_carcass <- DESeqDataSetFromMatrix(
  countData = count_matrix_TE_hypo_carcass, 
  colData = metadata_hypo_carcass,
  design = ~ Sex * Age..days.
)



##scale age as continuous ----
metadata_hypo_carcass_KW$Age..days. <- as.numeric(as.character(metadata_hypo_carcass_KW$Age..days.))

TE_dds_hypo_carcass_KW <- DESeqDataSetFromMatrix(
  countData = te_counts_filtered, 
  colData = metadata_hypo_carcass_KW,
  design = ~ Sex * Age..days.
)


TE_dds_hypo_carcass_KW <- DESeq(TE_dds_hypo_carcass_KW)

TE_vst <- varianceStabilizingTransformation(TE_dds_hypo_carcass_KW)


##results ----
TE_res_KW <- results(TE_dds_hypo_carcass_KW, name="Age..days.")

increasing_TE <- subset(TE_res_KW, 
                        padj < 0.05 & log2FoldChange > 0)

inc_TE_ns <- subset(TE_res_KW,
                    log2FoldChange > 0)

decreasing_TE <- subset(TE_res_KW, padj < 0.05 & log2FoldChange < 0)

dec_TE_ns <- subset(TE_res_KW,
                    log2FoldChange < 0)

TE_sex_res <- results(TE_dds_hypo_carcass_KW, name = "Sex_M_vs_F")

MvF_TE <- subset(TE_sex_res, 
                 padj <0.05)

