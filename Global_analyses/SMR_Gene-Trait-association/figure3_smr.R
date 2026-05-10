library(vroom)
library(tidyverse)
library(ComplexUpset)
library(ggplot2)


# 1/ Import gene-trait association SMR results (Supplementary Table 8)
smr_res <- vroom::vroom("ST8_SMR")

# Significant Sleep-trait SMR hits
signif_smr_sleep <- smr_res %>% 
  filter(!Trait %in% c("GGE", "Focal Epilepsy")) %>% 
  filter(p_SMR <= 0.05/10534)

# Significant Epilepsy SMR hits
signif_smr_epilepsy <- smr_res %>% 
  filter(Trait %in% c("GGE", "Focal Epilepsy")) %>% 
  filter(p_SMR <= 0.05/10534)

# 2/ Import Genes4Epilepsy (March 2026 version. Source https://bahlolab.github.io/Genes4Epilepsy/)
gene4epi <- vroom::vroom("Genes4Epilepsy_March2026.tsv")


# 3/ Import known clock genes
clockgenes <- c("CLOCK", "ARNTL1", "ARNTL2", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "RORA", "RORB", "RORC", "NR1D1", "NR1D2", "CSNK1D", "CSNK1E", "TIMELESS", "DBP", "TEF", "HLF", "BHLHE40", "BHLHE41")


# 4/ Prepare data for upset plot
# Create gene set for overlapping
genes <- unique(c(signif_smr_sleep$Gene,
                  signif_smr_epilepsy$Gene,
                  gene4epi$Gene),
                  clockgenes)

upset_df <- data.frame(
  Gene = genes,
  Sleep_SMR = genes %in% signif_smr_sleep$Gene,
  Epilepsy_SMR = genes %in% signif_smr_epilepsy$Gene,
  Genes4Epi = genes %in% gene4epi$Gene,
  Clockgene = genes %in% clockgenes) %>% 
  
  rename(
    `Sleep SMR hits` = "Sleep_SMR",
    `Epilepsy SMR hits` = "Epilepsy_SMR",
    `Genes4Epilepsy` = "Genes4Epi",
    `Clock genes` = "Clockgene")



# 5/ Find overlap genes for annotation
# Between Sleep_SMR and Epilepsy_SMR
overlap1 <- paste0(
  strrep("\n", 0),
  paste(intersect(signif_smr_sleep$Gene, signif_smr_epilepsy$Gene), collapse = "\n"))

# Between Gene4epi and Epilepsy_SMR
overlap2 <- intersect(gene4epi$Gene, signif_smr_epilepsy$Gene) %>% 
  paste(collapse = "\n") %>% 
  paste0(strrep("\n", 0))

# Between Gene4epi and Sleep_SMR
overlap3 <- intersect(gene4epi$Gene, signif_smr_sleep$Gene) %>% 
  paste(collapse = "\n") %>% 
  paste0(strrep("\n", 4))

# Between Gene4 epi and Clock genes
overlap4 <- paste0(
  strrep("\n", 0),
  paste(intersect(gene4epi$Gene, clockgenes), collapse = "\n"))



# 6/ Generate and Annotate Upset Plot
fig3_smr <- upset(
  upset_df,
  intersect = c("Sleep SMR hits", "Epilepsy SMR hits", "Genes4Epilepsy", "Clock genes"),
  base_annotations = 
    list(
      "Intersected Genes" = 
        intersection_size(counts = TRUE,
                          fill = "azure3",
                          color = "black",
                          text = list(size = 5),
                          text_colors = c(
                            on_bar = "white",
                            on_background = "black"),
                          bar_number_threshold = 1) +
        
        annotate("text", x = 4, y = 0.5, label = overlap1, size = 5, fontface = "italic") +
        annotate("text", x = 3, y = -0.6, label = overlap2, size = 5, fontface = "italic") +
        annotate("text", x = 2, y = 1.5, label = overlap4, size = 5, fontface = "italic") +
        annotate("text", x = 1, y = 4.5, label = overlap3, size = 5, fontface = "italic")),
  
  set_sizes = FALSE,
  min_degree = 2,
  matrix = intersection_matrix(
    geom = geom_point(size = 7),
    segment = geom_segment(size = 0.8))) +
  
  theme(
    axis.text.y = element_text(size = 15, color = "black"),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank())


ggsave(filename = file.path(plot_dir, "Fig3_SMR_Upset_plot.png"),
       plot = fig3_smr,
       width = 10,
       height = 8,
       units = "in",
       dpi = 999)
