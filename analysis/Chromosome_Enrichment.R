library(ggplot2)

main <- #dds objct
tester <- #gene_list

plot_name <- "Hypo_MF_Inc"

all_muller <- gsub("^[^-]+-([^-]+)-.*", "\\1", rownames(main))
deg_muller <- gsub("^[^-]+-([^-]+)-.*", "\\1", tester)

muller_levels <- sort(unique(all_muller))

ratio_df <- data.frame(
  muller    = muller_levels,
  bg_count  = as.numeric(table(all_muller)[muller_levels]),
  deg_count = as.numeric(table(factor(deg_muller, levels = muller_levels))[muller_levels])
)

ratio_df <- ratio_df[grepl("^Muller", ratio_df$muller), ]
ratio_df$ratio <- ratio_df$deg_count / ratio_df$bg_count  # moved to after filter

ggplot(ratio_df, aes(x = reorder(muller, -ratio), y = ratio)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = mean(ratio_df$ratio), linetype = "dashed", color = "grey40") +
  theme_classic() +
  labs(x = "Muller Element", y = "Proportion of genes that are Increasing")

ggsave(paste0("Results/Plots/Chromosome/", plot_name, ".pdf"))