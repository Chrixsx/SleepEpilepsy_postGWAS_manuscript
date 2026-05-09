library(vroom)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(ggrepel)
library(patchwork)


lava_bivar <- vroom::vroom("LAVA_bivariates")


# Import  2495 LAVA-defined loci (original source: https://github.com/josefin-werme/LAVA)
lava_loci <- vroom("Locus_blocks_s2500_m25_f1_w200.locfile") %>% 
  transmute(
    per_locus.id = LOC,
    chr  = as.integer(CHR),
    start = START,
    stop  = STOP
  ) %>%
  arrange(chr, start) %>%
  mutate(x_pos = row_number())


# Merge lava loci (all chromosme) with morningness-gge lava results
plot_df <- lava_loci %>% 
  left_join(lava_bivar, by = c("per_locus.id", "chr", "start", "stop"))


text_size = 35

p1 <- ggplot(plot_df,
               aes(x = x_pos, y = rho)) +
  geom_segment(aes(xend = x_pos, y = 0, yend = 0), colour = "grey90") +
  
  geom_point(aes(
      colour = ifelse(sig_annotate == "not_sig", "p_BH", `Trait 1`),
      shape  = `Trait 2`,
      alpha  = sig_annotate,
    ),
    size = 15, na.rm = TRUE
  ) +
  geom_hline(yintercept = 0, colour = "black", size = 3, linetype = "dashed") +
  scale_x_continuous(
    breaks = tapply(plot_df$x_pos, plot_df$chr, mean),
    labels = sort(unique(plot_df$chr))
  ) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.25)) +
  scale_color_manual(
    values = c(
      "Daytime Napping"    = "#9B58A5",
      "Daytime Sleepiness" = "#009E73",
      "Insomnia"           = "#D55E00",
      "Morningness"        = "#E69F00",
      "Sleep Duration"     = "#56B4E9",
      "Short Sleep"         = "#0072B2",
      "Long Sleep"        = "#B8D3DC"),
    
    breaks = c("Daytime Napping", "Daytime Sleepiness", "Insomnia", "Morningness", "Sleep Duration", "Short Sleep", "Long Sleep"),
    labels = c("Daytime Napping", "Daytime Sleepiness", "Insomnia", "Morningness", "Sleep Duration", "Short Sleep", "Long Sleep"),
    name = "Sleep-related traits (color): ") +
  
  scale_shape_manual(
    values = c("GGE" = 18, "Focal Epilepsy" = 20),
    name = "Epilepsy phenotypes (shape): ",
    na.value = NA,
    na.translate = FALSE) +
  
  scale_alpha_manual(
    values = c(
      "not_sig" = 0.4),  # grey, faint
    labels = c("not_sig" = "Not significant"),
    name = "Others: ",
    na.value = NA,
    na.translate = FALSE) +
  
  labs(x = "Chromosome", y = bquote("Correlation coefficient (local " * r[g] * ")")) +

  guides(
    shape  = guide_legend(order = 1, nrow = 1, byrow=F, override.aes = list(size = 10)),  # Epilepsy phenotype (1rd row)
    colour = guide_legend(order = 2, nrow =1, byrow=F, override.aes = list(size = 10)),  # Sleep trait (2nd row)
    alpha = "none") +
  
  theme_bw() +
  
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "horizontal",
    legend.box.just = "left",
    legend.title = element_text(size = 36, face = "bold"),
    legend.text  = element_text(size = 36),
    legend.title.align = 0,        
    legend.text.align  = 0,       
    
    legend.key.width   = unit(1, "lines"),
    legend.spacing.y  = unit(0.1, "lines"),
    
    legend.margin     = margin(t = 0.5, b = 0.5, l=0.5, r =0.5, unit = "cm"),
    panel.border      = element_rect(colour = "black", linewidth = 0.5),
    panel.background  = element_rect(fill = NA, colour = NA),
    plot.background   = element_rect(fill = NA, colour = NA),
    
    legend.box.background = element_rect(
      colour = "black",
      fill   = NA,
      linewidth = 0.2),
    
    legend.background = element_rect(fill = NA, colour = NA),
    legend.key        = element_rect(fill = NA, colour = NA),
    
    
    axis.text.x = element_text(color = "black", hjust = 1, size=text_size, angle = 0),
    axis.title.x = element_text(color = "black", size=text_size),
    
    axis.text.y = element_text(color = "black", size = text_size),
    axis.title.y = element_text(color = "black", size=text_size))


# Annotating the significant loci
p2 <- p1 +
  geom_text_repel(
    data = subset(lava_bivar_plot_df, sig_annotate == "p_BH" & !is.na(loci_band)),
    aes(label = loci_band),
    size = 10,
    box.padding   = 0.5,
    nudge_x       = +0.05,  # move label horizontally
    # nudge_y       = +0.2,  # move label vertically
    max.overlaps  = Inf,
    segment.color = "grey30",
    segment.size  = 0.3,
    show.legend   = FALSE)


# Add +ve and -ve correlation annotation
fig4b_miami <- p2 +
annotate("text", x = -1, y = 0.6, label = "+ve correlation\n\n\n",
         angle = 90, size = 12, fontface = "bold") +
annotate("text", x = -1, y = -0.6, label = "-ve correlation\n\n\n",
          angle = 90, size = 12, fontface = "bold")


# Export Image
ggsave(
  filename = "Figure4B_Miami_LAVA.png",
  plot     = fig4b_miami,
  width    = 35, height = 12,
  dpi      = 999)
