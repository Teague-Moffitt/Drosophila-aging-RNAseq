library(ggplot2)
library(DESeq2)
library(tidyverse)

# List one dds object per species
dds_list <- list(
  siamana    = dds_sia_carcass,
  hypocausta = dds_hypo_carcass_KW  
)

ages <- c("5", "10", "20", "40", "60")

plot_name <- "anti_color_sia_vs_hyp_shared"

genes <- c("maker-MullerB-snap-gene-44.37",
           "maker-MullerB-augustus-gene-204.102")

gene_names <- c("Ddc",
                "Catsup")

name_map <- setNames(gene_names, genes)

# extract long-format normalized counts + metadata for one dds
get_long_counts <- function(dds_obj, species_label) {
  normalized_counts <- counts(dds_obj, normalized = TRUE)
  counts_subset <- normalized_counts[rownames(normalized_counts) %in% genes, , drop = FALSE]
  
  counts_df <- as.data.frame(counts_subset)
  counts_df$gene <- name_map[rownames(counts_df)]
  
  counts_long <- counts_df %>%
    pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
    left_join(as.data.frame(colData(dds_obj)) %>%
                rownames_to_column("sample"),
              by = "sample") %>%
    mutate(Species = species_label)
  
  counts_long
}

# Combine all species into one long dataframe
counts_long_all <- purrr::imap_dfr(dds_list, get_long_counts)

# Summarize per gene / species / sex / age
counts_summary <- counts_long_all %>%
  group_by(gene, Species, Sex, Age..days.) %>%
  summarize(mean_expr = mean(expression),
            se = sd(expression) / sqrt(n()),
            .groups = "drop")

counts_summary$Age..days. <- factor(counts_summary$Age..days., levels = ages)

# --- Plot: color by Sex, linetype by Species, facet by gene ---
ggplot(counts_summary, aes(x = Age..days., y = mean_expr,
                           color = Sex, linetype = Species,
                           group = interaction(Sex, Species))) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_expr - se,
                    ymax = mean_expr + se), width = 0.2) +
  facet_wrap(~ gene, scales = "free_y") +
  theme_bw() +
  labs(x = "Age (days)", y = "Normalized Expression")

ggsave(paste0("Results/Plots/Single_Trajectories/", plot_name, ".pdf"))