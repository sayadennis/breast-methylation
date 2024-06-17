#' Dot plot to show most enriched results
#'
#' The input data frame should have an "estimate" and
#' a "FDR" columns.
#'
#' Top CG groups are determined by estimate (descending order).
#'
#' @param df data frame with the necessary columns
#' @param y the column to be plotted on y-axis
#' @param n number of hits to plot
#' @param order_by the column by which results are ordered
#' @param size_by the column by which dots are sized
#' @param color_by the column by which dots are colored
#' @param label_by the column for label (name of hit)
#' @param title plot title
#' @return grid plot object (by ggplot)
#'
#' @import ggplot2
#' @examples
#' dotPlot(data.frame(
#'   estimate = runif(10, 0, 10), FDR = runif(10, 0, 1), nD = runif(10, 10, 20),
#'   overlap = as.integer(runif(10, 0, 30)), group = "g", dbname = seq_len(10)
#' ))
#' @export
dotPlot <- function(df, y = "-log10(p)",
                    n = 10, order_by = "-log10(p)", title = "Pathway Enrichment Results",
                    label_by = "pathway", size_by = "overlap",
                    color_by = "enrichment", y_limits = NULL) {
  # Calculate -log10(p) if this is the specified y
  if (y == "-log10(p)") {
    df[["-log10(p)"]] <- -log10(df$p.value)
  }

  # Order dataframe by specified column
  df1 <- df[order(df[[order_by]]), ]
  df1 <- tail(df1, n = n) # tail because order is ascending

  # Ensure label_by is a factor with levels in the order of appearance
  df1[[label_by]] <- factor(df1[[label_by]], levels = df1[[label_by]])

  # Base plot
  p <- ggplot(df1) +
    geom_point(aes(.data[[label_by]], .data[[y]],
      size = .data[[size_by]], color = .data[[color_by]]
    )) +
    coord_flip() +
    ggtitle(title) +
    scale_color_gradient(low = "blue", high = "red") +
    ylab(y) +
    xlab("") +
    theme(axis.text.y = element_text(size = 12))

  if (!is.null(y_limits)) {
    p <- p + ylim(y_limits)
  }

  return(p)
}
