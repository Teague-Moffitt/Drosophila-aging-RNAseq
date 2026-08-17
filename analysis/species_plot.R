library(ggplot2)
library(DESeq2)
library(tidyr)
library(ggpubr)

plot_name <- "hyp2sia_M_aging_young_dotplot" #file name

res_sia <- results(
  dds_sia_carcass, 
  contrast = list(c("Age..days.","SexM.Age..days."))
) #result for species 1

res_hyp <- results(
  dds_hypo_carcass_KW_young,
  contrast = list(c("Age_c", "SexM.Age_c")) 
) # result for species 2

df_sia <- data.frame(gene = row.names(res_sia), sia = res_sia$log2FoldChange)
df_hyp <- data.frame(gene = row.names(res_hyp), hyp = res_hyp$log2FoldChange)

cross <- merge(df_sia, df_hyp, by = "gene") |> drop_na()

cross <- cross |> drop_na()

cor_test <- cor.test(cross$sia, cross$hyp, method = "spearman")

ggplot(
  data = cross,
  mapping = aes(x = sia, y = hyp)
  ) +
  geom_point(size = 1.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "firebrick") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  theme_sub_strip() +
  #xlim(-15,15) +
  #ylim(-30,30) +
  labs(
    title = "Transcriptome Wide Aging LogFC ",
    x = "D. siamana M Aging LogFC",
    y = "D. hypocausta M Aging (5,10,40) LogFC"
  ) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") 

ggsave(paste0("Results/Plots/dotplots/", plot_name, ".pdf"))
