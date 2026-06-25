library(ggplot2)
library(dplyr)
library(DESeq2)

# Get normalized counts for monotonically increasing genes
vsd_counts <- varianceStabilizingTransformation(dds_hypo_carcass_KW, blind=FALSE)

vsd_counts <- assay(vsd_hypo_carcass_KW)

# Subset to monotonic genes of interest
mono_genes <- mediate_hc[,1]  # or decreasing_ordered
mono_counts <- vsd_counts[mono_genes, ]
mono_counts <- mono_counts[, rownames(metadata_hypo_carcass_KW)]

# Get metadata
plot_meta <- metadata_hypo_carcass_KW

# Average expression per gene across age and sex
avg_expression <- as.data.frame(t(mono_counts)) %>%
  mutate(
    Sample = rownames(.),
    Sex = plot_meta[rownames(.), "Sex"],
    Age = as.numeric(as.character(plot_meta[rownames(.), "Age..days."]))
  ) %>%
  tidyr::pivot_longer(cols = all_of(mono_genes),
                      names_to = "Gene",
                      values_to = "Expression") %>%
  group_by(Sex, Age) %>%
  summarise(Mean_expression = mean(Expression),
            SE = sd(Expression)/sqrt(n()),
            .groups="drop")

# Plot
ggplot(avg_expression, aes(x=Age, y=Mean_expression, color=Sex)) +
  geom_line(linewidth=1) +
  geom_point(size=2) +
  geom_ribbon(aes(ymin=Mean_expression - SE,
                  ymax=Mean_expression + SE,
                  fill=Sex), alpha=0.2, color = NA) +
  scale_color_manual(values=c("F"="red", "M"="blue")) +
  scale_fill_manual(values=c("F"="red", "M"="blue")) +
  labs(#title="Average trajectory of monotonic decreasing genes",
       x="Age (days)",
       y="Mean VST expression") +
  theme_classic()

#ggsave("Plots/Avg_KW_mono_dec.pdf")