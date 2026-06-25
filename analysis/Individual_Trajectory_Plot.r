normalized_counts <- counts(dds_hypo_carcass_KW, normalized = TRUE)

genes <- head(mediate_hc[, 1], n = 79)

gene_names <- head(mediate_hc[, 3], n = 79)
  
counts_subset <- normalized_counts[rownames(normalized_counts) %in% genes, ]

library(tidyverse)

counts_df <- as.data.frame(counts_subset)
name_map <- setNames(gene_names, genes)
counts_df$gene <- name_map[rownames(counts_df)]

counts_long <- counts_df %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expression")

counts_long <- counts_long %>%
  left_join(as.data.frame(colData(dds_hypo_carcass)) %>%
              rownames_to_column("sample"),
            by = "sample")

counts_summary <- counts_long %>%
  group_by(gene, Sex, Age..days.) %>%
  summarize(mean_expr = mean(expression),
            se = sd(expression)/sqrt(n()))

counts_summary$Age..days. <- factor(counts_summary$Age..days.,
                                    levels = c("5", "10", "40", "60"))

ggplot(counts_summary, aes(x = Age..days., y = mean_expr,
                           color = Sex, group = Sex)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_expr - se,
                    ymax = mean_expr + se), width = 0.2) +
  facet_wrap(~ gene, scales = "free_y") +
  theme_bw() +
  labs(x = "Age (days)", y = "Normalized Expression")