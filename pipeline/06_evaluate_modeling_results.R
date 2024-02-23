library(ggsankey)
library(ggplot2)
library(dplyr)

space <- 30000

din <- "/projects/p30791/methylation/sesame_out/differential_methylation"

trends <- read.csv(paste0(din, "/dml_up_down_pairwise_trends.csv"))

df <- trends %>% make_long("All.probes", "Normal.vs.CUB", "CUB.vs.OQ", "OQ.vs.AN", "AN.vs.TU")

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
  scale_fill_viridis_d() +
  theme_sankey(base_size = 16) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  ) +
  scale_x_discrete(labels = c("All probes", "Normal vs CUB", "CUB vs OQ", "OQ vs AN", "AN vs TU"))

plot_filepath <- "/projects/p30791/methylation/plots/dml_sankey.png"
ggsave(plot_filepath, plot = p, width = 10, height = 6)
