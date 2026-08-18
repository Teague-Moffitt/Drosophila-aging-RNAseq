library(DESeq2)
library(here)
#ALL DATA ----

##make metadata ----
metadata <- read.csv(
  #metadata.csv, 
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
count_dir <- #read count directory

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

#center sex at mean age

metadata_hypo_carcass_KW$Age_c <- metadata_hypo_carcass_KW$Age..days. - mean(metadata_hypo_carcass_KW$Age..days., na.rm = TRUE)


##subset matrix ----
count_matrix_hypo_carcass_KW <- count_matrix[, hypo_carcass_KW_samples]



##make DDS object hypocarcass ----
dds_hypo_carcass_KW <- DESeqDataSetFromMatrix(
  countData = count_matrix_hypo_carcass_KW, 
  colData = metadata_hypo_carcass_KW,
  design = ~ Sex * Age_c
)


dds_hypo_carcass_KW <- DESeq(dds_hypo_carcass_KW)

vsd_hypo_carcass_KW <- vst(dds_hypo_carcass_KW, blind = TRUE)


##Monotonic Genes ----
inc_hypo_f <- subset(results(dds_hypo_carcass_KW, name = "Age_c"), padj < 0.05 & log2FoldChange > 0)

dec_hypo_f <- subset(results(dds_hypo_carcass_KW, name = "Age_c"), padj < 0.05 & log2FoldChange < 0)

res_sex_hyp <- results(dds_hypo_carcass_KW, name = "Sex_M_vs_F")

res_age_female <- results(dds_hypo_carcass_KW, name = "Age_c")

res_age_male <- results(
  dds_hypo_carcass_KW,
  contrast = list(c("Age_c", "SexM.Age_c"))
)

inc_hypo_m <- subset(res_age_male, padj < 0.05 & log2FoldChange > 0)

dec_hypo_m <- subset(res_age_male, padj < 0.05 & log2FoldChange < 0)

#D. SIAMANA liftoff genome ----

##make metadata ----

metadata_sia <- subset(
  metadata, 
  Species == "D. siamana"
)

##create list of files ----
count_dir_sia <- here("Data", "Dsia_read_counts")

files <- list.files(
  path = count_dir_sia,
  pattern = "\\.txt$",
  full.names = TRUE
)


##match samples ----
sample_ids_sia <- sub(
  "_counts\\.txt$",
  "",
  basename(files)
)

metadata_order_sia <- rownames(metadata_sia)

files <- files[
  match(metadata_order_sia, sample_ids_sia)
]


##create count matrix ----
count_list_sia <- lapply(files, function(f) {
  
  df_sia <- read.table(
    f, 
    header = TRUE,
    comment.char = "#",
    stringsAsFactors = FALSE
  )
  
  counts_sia <- df_sia[, c(1, ncol(df_sia))]
  
  sample_name_sia <- sub(
    "_counts\\.txt$",
    "",
    basename(f)
  )
  
  colnames(counts_sia) [2] <- sample_name_sia
  
  counts_sia
})


##merge matrix ----
count_matrix_sia <- Reduce(function(x, y) {
  merge(x, y, by = "Geneid")
}, count_list_sia)

rownames(count_matrix_sia) <- count_matrix_sia$Geneid

count_matrix_sia <- count_matrix_sia[, -1]

count_matrix_sia <- as.matrix(count_matrix_sia)


##make dds object ----
dds_sia <- DESeqDataSetFromMatrix(
  countData = count_matrix_sia, 
  colData = metadata_sia,
  design = ~ Sex + Tissue + Age..days. 
)

dds_sia <- DESeq(dds_sia)
vsd_sia <- vst(dds_sia, blind = TRUE)

#D. SIAMANA CARCASS KW0 ----

##create sia carcass KW metadata ----
metadata_sia_carcass <- subset(
  metadata, 
  Tissue == "Carcass" & Species == "D. siamana"
)

sia_carcass_samples <- rownames(metadata_sia_carcass)

metadata_sia_carcass$Age..days. <- as.numeric(as.character(metadata_sia_carcass$Age..days.))



##subset matrix ----
count_matrix_sia_carcass <- count_matrix_sia[, sia_carcass_samples]



##make DDS object siacarcass ----
dds_sia_carcass <- DESeqDataSetFromMatrix(
  countData = count_matrix_sia_carcass, 
  colData = metadata_sia_carcass,
  design = ~ Sex * Age..days.
)


dds_sia_carcass <- DESeq(dds_sia_carcass)

vsd_sia_carcass <- vst(dds_sia_carcass, blind = TRUE)


##Monotonic Genes ----
inc_sia_f <- subset(results(dds_sia_carcass, name = "Age..days."), padj < 0.05 & log2FoldChange > 0)

dec_sia_f <- subset(results(dds_sia_carcass, name = "Age..days."), padj < 0.05 & log2FoldChange < 0)

res_sia_age_female <- results(dds_sia_carcass, name = "Age..days.")

res_sia_age_male <- results(
  dds_sia_carcass,
  contrast = list(c("Age..days.", "SexM.Age..days."))
)

inc_sia_m <- subset(res_sia_age_male, padj < 0.05 & log2FoldChange > 0)

dec_sia_m <- subset(res_sia_age_male, padj < 0.05 & log2FoldChange < 0)

#D. hypocausta germline ----

##create hypo gl metadata ----
metadata_hypo_gl <- subset(
  metadata, 
  Tissue %in% c("Ovaries", "Testes") & Species == "D. hypocausta"
)

hypo_gl_samples <- rownames(
  metadata_hypo_gl[
    grepl(
      "^KW",
      rownames(metadata_hypo_gl)
    ), 
  ]
)

metadata_hypo_gl <- metadata_hypo_gl[hypo_gl_samples, ]

metadata_hypo_gl$Age..days. <- as.numeric(as.character(metadata_hypo_gl$Age..days.))

#center sex at mean age

metadata_hypo_gl$Age_c <- metadata_hypo_gl$Age..days. - mean(metadata_hypo_gl$Age..days., na.rm = TRUE)


##subset matrix ----
count_matrix_hypo_gl <- count_matrix[, hypo_gl_samples]



##make DDS object hypo germline ----
dds_hypo_gl <- DESeqDataSetFromMatrix(
  countData = count_matrix_hypo_gl, 
  colData = metadata_hypo_gl,
  design = ~ Sex * Age_c
)


dds_hypo_gl <- DESeq(dds_hypo_gl)

vsd_hypo_gl <- vst(dds_hypo_gl, blind = TRUE)


##Monotonic Genes ----
inc_hypo_gl_f <- subset(results(dds_hypo_gl, name = "Age_c"), padj < 0.05 & log2FoldChange > 0)

dec_hypo_gl_f <- subset(results(dds_hypo_gl, name = "Age_c"), padj < 0.05 & log2FoldChange < 0)

res_sex_hypo_gl <- results(dds_hypo_gl, name = "Sex_M_vs_F")

res_hypo_gl_age_female <- results(dds_hypo_gl, name = "Age_c")

res_hypo_gl_age_male <- results(
  dds_hypo_gl,
  contrast = list(c("Age_c", "SexM.Age_c"))
)

inc_hypo_gl_m <- subset(res_hypo_gl_age_male, padj < 0.05 & log2FoldChange > 0)

dec_hypo_gl_m <- subset(res_hypo_gl_age_male, padj < 0.05 & log2FoldChange < 0)

#D. siamana germline ----

##create sia gl metadata ----
metadata_sia_gl <- subset(
  metadata, 
  Tissue %in% c("Ovaries", "Testes") & Species == "D. siamana"
)

sia_gl_samples <- rownames(metadata_sia_gl)

metadata_sia_gl <- metadata_sia_gl[sia_gl_samples, ]

metadata_sia_gl$Age..days. <- as.numeric(as.character(metadata_sia_gl$Age..days.))

#center sex at mean age

metadata_sia_gl$Age_c <- metadata_sia_gl$Age..days. - mean(metadata_sia_gl$Age..days., na.rm = TRUE)


##subset matrix ----
count_matrix_sia_gl <- count_matrix[, sia_gl_samples]



##make DDS object sia germline ----
dds_sia_gl <- DESeqDataSetFromMatrix(
  countData = count_matrix_sia_gl, 
  colData = metadata_sia_gl,
  design = ~ Sex * Age_c
)


dds_sia_gl <- DESeq(dds_sia_gl)

vsd_sia_carcass_KW <- vst(dds_sia_gl, blind = TRUE)


##Monotonic Genes ----
inc_sia_gl_f <- subset(results(dds_sia_gl, name = "Age_c"), padj < 0.05 & log2FoldChange > 0)

dec_sia_gl_f <- subset(results(dds_sia_gl, name = "Age_c"), padj < 0.05 & log2FoldChange < 0)

res_sex_sia_gl <- results(dds_sia_gl, name = "Sex_M_vs_F")

res_sia_gl_age_female <- results(dds_sia_gl, name = "Age_c")

res_sia_gl_age_male <- results(
  dds_sia_gl,
  contrast = list(c("Age_c", "SexM.Age_c"))
)

inc_sia_gl_m <- subset(res_sia_gl_age_male, padj < 0.05 & log2FoldChange > 0)

dec_sia_gl_m <- subset(res_sia_gl_age_male, padj < 0.05 & log2FoldChange < 0)

#D. HYPOCASTA CARCASS KW w/out 60 day old flies ----

##create hypo carcass KW young metadata ----
metadata_hypo_carcass_young <- subset(
  metadata, 
  Tissue == "Carcass" & Species == "D. hypocausta" & Age..days. < 60
)

hypo_carcass_KW_samples_young <- rownames(
  metadata_hypo_carcass_young[
    grepl(
      "^KW",
      rownames(metadata_hypo_carcass_young)
    ), 
  ]
)

metadata_hypo_carcass_KW_young <- metadata_hypo_carcass_young[hypo_carcass_KW_samples_young, ]

metadata_hypo_carcass_KW_young$Age..days. <- as.numeric(as.character(metadata_hypo_carcass_KW_young$Age..days.))

#center sex at mean age

metadata_hypo_carcass_KW_young$Age_c <- metadata_hypo_carcass_KW_young$Age..days. - mean(metadata_hypo_carcass_KW_young$Age..days., na.rm = TRUE)


##subset matrix ----

count_matrix_hypo_carcass_KW_young <- count_matrix[, hypo_carcass_KW_samples_young]



##make DDS object hypocarcass ----
dds_hypo_carcass_KW_young <- DESeqDataSetFromMatrix(
  countData = count_matrix_hypo_carcass_KW_young, 
  colData = metadata_hypo_carcass_KW_young,
  design = ~ Sex * Age_c
)


dds_hypo_carcass_KW_young <- DESeq(dds_hypo_carcass_KW_young)


#Repeat Expression ----

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

count_matrix_TE_sia_carcass <- count_matrix_TE[, sia_carcass_samples]

##make DEseq object ----
### for all ----

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

TE_dds_sia_carcass <- DESeqDataSetFromMatrix(
  countData = count_matrix_TE_sia_carcass, 
  colData = metadata_sia_carcass,
  design = ~ Sex * Age..days.
)

##scale age as continuous ----
metadata_hypo_gl$Age..days. <- as.numeric(as.character(metadata_hypo_gl$Age..days.))

TE_dds_hypo_carcass_KW <- DESeqDataSetFromMatrix(
  countData = te_counts_filtered, 
  colData = metadata_hypo_carcass_KW,
  design = ~ Sex * Age..days.
)


TE_dds_hypo_carcass_KW <- DESeq(TE_dds_hypo_carcass_KW)

TE_vst <- varianceStabilizingTransformation(TE_dds_hypo_carcass_KW)


##results ----
TE_res_KW <- results(TE_dds_hypo_carcass_KW, name="Age..days.")

inc_TE_f <- subset(TE_res_KW, 
                   padj < 0.05 & log2FoldChange > 0)

inc_TE_ns <- subset(TE_res_KW,
                    log2FoldChange > 0)

dec_TE_f <- subset(TE_res_KW, padj < 0.05 & log2FoldChange < 0)

dec_TE_ns <- subset(TE_res_KW,
                    log2FoldChange < 0)

TE_sex_res <- results(TE_dds_hypo_carcass_KW, name = "Sex_M_vs_F")

MvF_TE <- subset(TE_sex_res, 
                 padj <0.05)

TE_res_male <- results(
  TE_dds_hypo_carcass_KW,
  contrast = list(c("Age..days.", "SexM.Age..days."))
)

inc_TE_m <- subset(TE_res_male, 
                   padj < 0.05 & log2FoldChange > 0)

dec_TE_m <- subset(TE_res_male,
                   padj < 0.05 & log2FoldChange < 0)
