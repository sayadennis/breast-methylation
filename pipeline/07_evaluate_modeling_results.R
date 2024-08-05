library(ggsankey)
library(ggplot2)
library(dplyr)

space <- 30000

din <- "/projects/p30791/methylation/differential_methylation"
plot_dir <- "/projects/p30791/methylation/plots"

trends <- read.csv(paste0(din, "/dml_hyper_hypo_pairwise_trends.csv"))

df <- trends %>% make_long("All.probes", "CFN.vs.CUB", "CUB.vs.OQ", "OQ.vs.AN", "AN.vs.TU")

node_colors <- c(
  "hyper" = "#FF00FF", # Magenta
  "hypo" = "#00FFFF", # Cyan
  "n.d." = "#808080" # Gray
)

p <- ggplot(df, aes(
  x = x,
  next_x = next_x,
  node = node,
  next_node = next_node,
  label = node,
  fill = factor(node)
)) +
  geom_sankey(flow.alpha = 0.75, node.color = 1, flow.color = 1, flow.fill = 8, space = space) +
  geom_sankey_label(size = 3.5, color = 1, fill = "white", space = space) +
  scale_fill_manual(values = node_colors) + # Use custom colors
  theme_sankey(base_size = 16) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  ) +
  scale_x_discrete(labels = c("All probes", "CFN vs CUB", "CUB vs OQ", "OQ vs AN", "AN vs TU"))

plot_filepath <- paste0(plot_dir, "/dml_sankey.png")
ggsave(plot_filepath, plot = p, width = 7, height = 4.2)
