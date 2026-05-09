library(ggplot2)

fig2_ldsc <- ggcorrplot(
  corr.matrix,
  type = "full",
  method = "square",
  hc.order = FALSE,
  legend.title = expression("global" ~ r[g]),
  lab = TRUE,
  lab_col = "black",
  lab_size = 15,
  outline.color = "black",
  colors = c("blue", "white", "red"),
  ggtheme = theme_bw() +
    theme(
      axis.text.x = element_text(color="black", size=20),
      axis.text.y = element_text(color = "black"),
      
      legend.text = element_text(size = 25),
      legend.key.size = unit(1.5, "cm"),
      legend.title = element_text(size = 30),
      legend.title.align = .5,
      legend.position = "right", # keep figure legend to right hand side
      legend.justification = c(1, 0.5),
      legend.spacing.y = unit(0, "pt"), # adjust spacing in figure legend
      legend.margin = margin(t = 10, r = 10, b = 10, l = 10),
      
      plot.margin = margin(t=2, b=2, l=0,r=0, unit="cm"))) +
  
  guides(fill = guide_colorbar(barwidth = 3.0, barheight = 17))


# Annotating significant results
annotated_fig2_ldsc <- fig2_ldsc +
  annotate("text", x=4, y=1.30, label= "✶", col="red4", size=20) +
  annotate("text", x=2, y=1.30, label= "✶", col="red4", size=20)


# Export image
ggsave(filename = "imagestorage_path/figure2_LDSC_glocal_correlation.png",
       plot = annotated_fig2_ldsc,
       width = 50,
       height = 30,
       units = "cm",
       dpi=999)
