library(ggplot2)
## PCA plot of all data: species, tissue ----
plotPCA(TE_vsd, intgroup = c("Tissue", "Species")) + ggtitle("D. hypocausta and D. Siamana Differential Expression")

ggsave("Plots/all_tissue_species.pdf")

## PCA plot of hypo data: tissue ----
plotPCA(vsd_hypo, intgroup = c("Tissue", "Sex")) + ggtitle("D. hypocausta Ovarie, Teste and Carcass Differential Expression")

ggsave("Plots/hypo_tissue.pdf")

##PCA plot of hypo_carcass Sex vs Age (days) Interaction ----

plotPCA(TE_vst, intgroup = c("Sex", "Age..days."))

pcaData <- plotPCA(TE_vst,
                   intgroup = c("Sex", "Age..days."),
                   returnData = TRUE)

pcaData$Age..days. <- factor(
  pcaData$Age..days.,
  levels = sort(as.numeric(as.character(unique(pcaData$Age..days.))))
)

library(ggplot2)

ggplot(
  pcaData,
  aes(PC1, PC2,
      shape = Sex,
      color = Age..days.)
) +
  geom_point(size = 4, alpha = 0.9) +
  labs(
    # title = "D. hypocausta Carcass Differential Expression",
    shape = "Sex",
    color = "Age (days)"
  ) +
  scale_color_viridis_d()
+ scale_shape_discrete()
 ggsave("Plots/Dhypo_carcass_sex_age.pdf")

## PCA plot of all TE data: species, tissue ----
plotPCA(TE_vsd, intgroup = c("Tissue", "Species")) + ggtitle("D. hypocausta and D. Siamana Differential TE Expression")

ggsave("Plots/TE_all_tissue_species.pdf")

##PCA plot of TE_hypo data: tissue ----
plotPCA(TE_vsd_hypo, intgroup = c("Tissue", "Sex")) + ggtitle("D. hypocausta Ovarie, Teste and Carcass Differential TE Expression")

ggsave("Plots/TE_hypo_tissue.pdf")

