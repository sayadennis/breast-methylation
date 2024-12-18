library(ggplot2)

source("/home/srd6051/breast-methylation/exploration/dotPlot.R")

din <- "/projects/p30791/methylation/differential_methylation"

gometh_results <- list(
  "TFBS_HDB_vs_CUB_only" = read.csv(
    paste0(din, "/gometh_KEGG_TFBS_HDB_vs_CUB_only.csv")
  ),
  "TFBS_AN_vs_TU_only" = read.csv(
    paste0(din, "/gometh_KEGG_TFBS_AN_vs_TU_only.csv")
  ),
  "TFBS_HDBvCUB_ANvTU_overlap" = read.csv(
    paste0(din, "/gometh_KEGG_TFBS_HDBvCUB_ANvTU_overlap.csv")
  )
)

pathway_sizes <- c(
  gometh_results[["TFBS_HDB_vs_CUB_only"]]$N,
  gometh_results[["TFBS_AN_vs_TU_only"]]$N
)
overlaps <- c(
  gometh_results[["TFBS_HDB_vs_CUB_only"]]$DE,
  gometh_results[["TFBS_AN_vs_TU_only"]]$DE
)

longname <- "Inflammatory mediator regulation of TRP channels"
shortname <- "Inflammatory mediator regulation"
gometh_results[["TFBS_HDB_vs_CUB_only"]][
  gometh_results[["TFBS_HDB_vs_CUB_only"]]$Description == longname, "Description"
] <- shortname

n_show <- 7

for (setname in names(gometh_results)) {
  test <- gometh_results[[setname]]
  # Plot dot plot
  colnames(test) <- c("HSA", "Pathway", "Pathway size", "Overlap", "p.value", "q.value")
  p <- dotPlot(
    test,
    n = n_show,
    color_by = "Pathway size", y_limits = c(1.25, 3.0),
    label_by = "Pathway", size_by = "Overlap", title = "",
  )
  ggsave(
    filename = paste0("/home/srd6051/dotplot_gometh_", setname, ".png"),
    plot = p,
    width = 24,
    height = 2 + 0.8 * n_show,
    units = "cm"
  )
}
