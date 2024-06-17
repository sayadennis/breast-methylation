library(ggplot2)

dotPlot <- function(df, y = "-log10(p)",
                    n = 10, order_by = "-log10(p)", title = "Pathway Enrichment Results",
                    label_by = "pathway", size_by = "overlap", color_by = "enrichment") {
  # Calculate -log10(p) if this is the specified y
  if (y == "-log10(p)") {
    df[["-log10(p)"]] <- -log10(df$p.value)
  }

  # Order dataframe by specified column
  df1 <- df[order(df[[order_by]]), ]
  df1 <- head(df1, n = n)

  # Ensure label_by is a factor with levels in the order of appearance
  df1[[label_by]] <- factor(df1[[label_by]], levels = df1[[label_by]])

  # Return plot
  ggplot(df1) +
    geom_point(aes(.data[[label_by]], .data[[y]],
      size = .data[[size_by]], color = .data[[color_by]]
    )) +
    coord_flip() +
    ggtitle(title) +
    scale_color_gradient(low = "blue", high = "red") +
    ylab(y) +
    xlab("")
}

## Create dummy data to test plotting

# Set seed for reproducibility
set.seed(123)

# Generate dummy data
Pathway_ID <- 1:12
Pathway_name <- paste("Pathway", LETTERS[1:12])
p_value <- round(runif(12, min = 0.0001, max = 0.05), 5)
q_value <- round(runif(12, min = 0.0001, max = 0.05), 5)
overlap <- sample(5:50, 12, replace = TRUE)
enrichment <- round(runif(12, min = 1, max = 2.5), 2)

# Create data frame
dummy_data <- data.frame(
  `Pathway ID` = Pathway_ID,
  `Pathway name` = Pathway_name,
  `p.value` = p_value,
  `q.value` = q_value,
  `overlap` = overlap,
  `enrichment` = enrichment
)

print(dummy_data)

## Try plotting

p <- dotPlot(
  dummy_data,
  y = "-log10(p)",
  n = 10, # number of hits to show
  title = "Pathway Enrichment Results",
  order_by = "-log10(p)",
  color_by = "enrichment",
  label_by = "Pathway.name",
  size_by = "overlap"
)

plot_dir <- "/home/srd6051"

ggsave(
  filename = paste0(plot_dir, "/dummy_dotplot.png"),
  plot = p,
  width = 16,
  height = 2 + 0.8 * dim(dummy_data)[1],
  units = "cm"
)
