library(ggplot2)
library(dplyr)
library(DESeq2)

metadat_avg <- #metadata
gene_list <-#gene_names
given_dds <- #dds object

plot_name <- #file name
plot_title <- #title
  
# Get normalized counts 
vsd_counts <- varianceStabilizingTransformation(given_dds, blind=FALSE)

vsd_counts <- assay(vsd_counts)

# Subset to  genes of interest
mono_genes <- (gene_list)  
mono_genes <- intersect(mono_genes, rownames(vsd_counts))
mono_counts <- vsd_counts[mono_genes,]
mono_counts <- mono_counts[, rownames(metadat_avg)]
common_samples <- intersect(colnames(mono_counts), rownames(metadat_avg))
mono_counts <- mono_counts[, common_samples]

# Get metadata
plot_meta <- metadat_avg[common_samples, ]

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
  labs(title= plot_title,
       x="Age (days)",
       y="Mean VST expression") +
  theme_classic()

ggsave(paste0("Results/Plots/Avg_Trajectories/", plot_name, ".pdf"))