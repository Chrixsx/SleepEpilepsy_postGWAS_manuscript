library(vroom)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(igraph)
library(scales)

# Import bivariate LAVA correlation results
lava_bivar <- vroom::vroom("LAVA_bivariates") %>% 
  mutate(
    pair = paste0(`Trait 1`, "-", `Trait 2`)) %>% 
  group_by(pair) %>% 
  mutate(
    shared_loci = list(Cytogenic_band),
    weight = n())

# Create connection between significant sleep-epilepsy pair, and count number of correlation
edges_df <- lava_bivar %>%
  mutate(
    loci_label = map_chr(shared_loci, ~ paste(.x, collapse = ", "))
  ) %>%
  rename(sleep = `Trait 1`, epilepsy = `Trait 2`) %>% 
  select(sleep, epilepsy, weight, loci_label, pair)

sleep_traits <- unique(edges_df$sleep)
epi_traits   <- unique(edges_df$epilepsy)

# Prepare Node for network
nodes <- tibble(
  name  = c(sleep_traits, epi_traits),
  group = c(rep("Sleep-related traits", length(sleep_traits)),
            rep("Epilepsy phenotypes", length(epi_traits)))
) %>%
  mutate(type = group == "Epilepsy phenotypes")

# Prepare Edge for network
edges_graph <- edges_df %>%
  transmute(
    from       = sleep,
    to         = epilepsy,
    weight,
    loci_label)

# Plot undirected network plot from Nodes and Edges
g1 <- tbl_graph(nodes = nodes, edges = edges_graph, directed = FALSE)

lay <- layout_as_bipartite(g1, hgap = 0.5, vgap = 0.3)


layout_vertical <- data.frame(
  x = lay[,2] * 0.3,  # Reduce horizontal spread
  y = lay[,1] * 0.2)

# Raw plot
p1 <- ggraph(g1, layout = layout_vertical)

# Adjusting edge (weight) size and color
p2 <- p1 +
  geom_edge_link(
    aes(width = weight),
    colour = "grey45",
    alpha  = 0.7
  ) +
  scale_edge_width(
    range  = c(3, 10),
    name   = "Number of significantly correlated loci",
    labels = scales::label_number(accuracy = 1))


# Define data shape of epilepsy types and color for sleep-related traits
p3 <- p2 +
  geom_node_point(aes(color = ifelse(group == "Sleep-related traits", name, "black"), 
                      shape = ifelse(group == "Epilepsy phenotypes", name, 19)), 
                  size = 39) +
  # Removed geom_node_text layers for manual annotation
  scale_color_manual(values = c(
    "Daytime Napping"    = "#9B58A5",
    "Daytime Sleepiness" = "#009E73",
    "Insomnia"           = "#D55E00",
    "Morningness"        = "#E69F00",
    "Sleep Duration"     = "#56B4E9",
    "Short Sleep"        = "#0072B2",
    "Long Sleep"         = "#B8D3DC"
  ), na.value = "grey50") +
  scale_shape_manual(values = c("GGE" = 18, "Focal Epilepsy" = 20),
                     na.value = 19) +
  guides(color = "none", shape = "none") +
  theme_graph()


# Making figure legend box
p4 <-  p3 +
  coord_cartesian(clip = "off") +  # Prevents ggplot clipping
  theme_graph() +
  theme(
    legend.position = "bottom",                # Move to bottom
    legend.direction = "horizontal",           # Horizontal layout
    legend.margin = margin(t = 1, b = 1),      # Top/bottom padding
    legend.title = element_text(size = 45),
    legend.text = element_text(size = 45),
    legend.key.size = unit(2.5, "cm"),

    plot.margin = margin(t = 3,  b = 0.5, r = 19.5, l = 14.5, unit = "cm"))


# Annotate each data point
text_size = 27
x_gap = 0.103

fig4a_network <-
  p4 +
  annotate("text", x = x_gap, y = +0.6, label = "Long Sleep", size = text_size, color = "black", hjust = 0) +
  annotate("text", x = x_gap, y = +0.5, label = "Short Sleep", size = text_size, color = "black", hjust = 0) +
  annotate("text", x = x_gap, y = +0.4, label = "Morningness", size = text_size, color = "black", hjust = 0) +
  annotate("text", x = x_gap, y = +0.3, label = "Daytime Sleepiness", size = text_size, color = "black", hjust = 0) +
  annotate("text", x = x_gap, y = +0.2, label = "Daytime Napping", size = text_size, color = "black", hjust = 0) +
  annotate("text", x = x_gap, y = +0.1, label = "Sleep Duration", size = text_size, color = "black", hjust = 0) +
  annotate("text", x = x_gap, y = +0.0, label = "Insomnia", size = text_size, color = "black", hjust = 0) +
 
  annotate("text", x = -0.012, y = +0.55, label = "Focal Epilepsy", size = text_size, color = "black", hjust = 1) + 
  annotate("text", x = -0.012, y = +0.3, label = "GGE", size = text_size, color = "black", hjust = 1) 


# Export Image
ggsave("figure4A_LAVA_Network.png",
       plot = fig4a_network,
       width = 20, height = 18, dpi = 999)
