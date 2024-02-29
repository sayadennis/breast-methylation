library(stringr)
library(sesame)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Please provide 3 arguments: ",
    "1) input directory, ",
    "2) output directory to write CSV files to, and ",
    "3) output directory for plots.",
    call. = FALSE
  )
} else {
  din <- args[1] # e.g. "/projects/p30791/methylation/sesame_out/differential_methylation" # nolint
  dout <- args[2] # e.g. "/projects/p30791/methylation/sesame_out/KYCG"
  plot_dir <- args[3] # e.g. "/projects/p30791/methylation/plots/KYCG"
}

if (!file.exists(dout)) {
  dir.create(dout)
  cat("Directory created:", dout, "\n")
}

## Read the query probe sets
probeset_names <- c(
  "hypo_in_TU_HER2",
  "hypo_in_TU_ER",
  "hyper_in_TU_ER",
  "CUB_down_from_Normal",
  "CUB_up_from_Normal",
  "AN_up_from_OQ",
  "AN_down_from_OQ",
  "AN_up_and_TU_down",
  "AN_up_and_TU_nd_or_up",
  "AN_down_and_TU_nd_or_down",
  "AN_down_and_TU_up",
  "TU_down_from_AN",
  "TU_up_from_AN"
)

probe_sets <- list()

for (setname in probeset_names) {
  file_path <- paste0(din, "/probe_set_", setname, ".txt")
  probe_sets[[setname]] <- readLines(file_path)
}

## Download public Databases
dbs <- KYCG_listDBGroups("EPIC")
dbs <- dbs[str_starts(dbs$Title, fixed("KYCG.EPIC.")), ]
db_list <- list()
for (db_id in dbs$Title) {
  db_list[[db_id]] <- KYCG_getDBs(db_id)
}

########################################################
#### Categorical vs Categorical Enrichment Analysis ####
########################################################

dbs_to_plot <- c(
  "chromHMM",
  "HMconsensus",
  # "seqContext",
  "TFBSconsensus"
)

for (setname in probeset_names) {
  print(paste0("Testing ", setname, "..."))
  query <- probe_sets[[setname]]

  # Test enrichment for public databases
  for (db_id in dbs$Title) {
    results <- testEnrichment(query, db_list[[db_id]], platform = "EPIC")

    # Select statistically significant results
    db_name <- strsplit(db_id, "\\.")[[1]][3]
    results <- results[results$p.value < 0.05, ]
    if (dim(results)[1] == 0) {
      # only save results if any are significant
      next
    }

    # Fix rows where p-value and FDR-adjusted p-values suffer from underflow
    min_nonzero_value <- min(results$p.value[results$p.value > 0], na.rm = TRUE)
    results$p.value[results$p.value == 0] <- min_nonzero_value

    min_nonzero_value <- min(results$FDR[results$FDR > 0], na.rm = TRUE)
    results$FDR[results$FDR == 0] <- min_nonzero_value

    file_path <- paste0(dout, "/testEnrichment_", setname, "_", db_name, ".csv")
    write.csv(results, file_path, row.names = FALSE)
    if (db_name %in% dbs_to_plot) {
      dotplot <- KYCG_plotDot(
        results,
        n = 20,
        title = paste0(
          "Enriched ", gsub("consensus", "", db_name), " -- ",
          gsub("_", " ", setname)
        )
      )
      plot_filepath <- paste0(
        plot_dir, "/dotplot_", setname, "_", db_name, ".png"
      )
      num_sig_results <- dim(results)[1]
      if (num_sig_results >= 20) {
        height <- 5
      } else {
        height <- 3.5 + num_sig_results * 1.5 / 20
      }
      ggsave(plot_filepath, plot = dotplot, width = 6, height = height)
    }
  }

  # Test enrichment for genes
  results <- testEnrichment(
    query,
    KYCG_buildGeneDBs(query, max_distance = 100000, platform = "EPIC"),
    platform = "EPIC"
  )
  results <- results[results$p.value < 0.05, ]
  if (dim(results)[1] > 0) {
    # Fix rows where p-value and FDR-adjusted p-values suffer from underflow
    min_nonzero_value <- min(results$p.value[results$p.value > 0], na.rm = TRUE)
    results$p.value[results$p.value == 0] <- min_nonzero_value

    min_nonzero_value <- min(results$FDR[results$FDR > 0], na.rm = TRUE)
    results$FDR[results$FDR == 0] <- min_nonzero_value

    file_path <- paste0(dout, "/testEnrichment_", setname, "_genes.csv")
    write.csv(results, file_path, row.names = FALSE)
  }
}

########################################
#### GO/Pathway Enrichment Analysis ####
########################################

regs <- sesameData_getTxnGRanges("hg38")

source("sesameData_annoProbes_custom.R")

library(gprofiler2)

for (setname in probeset_names) {
  print(paste0("Working on ", setname, "..."))
  query <- probe_sets[[setname]]
  genes <- sesameData_annoProbes_custom(
    query, regs,
    platform = "EPIC", return_ov_features = TRUE
  )
  gostres <- gost(
    genes$gene_name,
    organism = "hsapiens"
  )
  if (!is.null(gostres)) {
    df <- gostres$result
    df <- df[, !(names(df) %in% "parents")]
    file_path <- paste0(dout, "/gost_result_", setname, ".csv")
    write.csv(df[order(df$p_value), ], file_path, row.names = FALSE)
  } else {
    print(paste0("No GO/Pathway Enrichment results for", setname))
  }
}

# nolint start: commented_code_linter.

# ###########################################################
# #### Categorical vs Continuous Set Enrichment Analysis ####
# ###########################################################
#
# # Query can be continuous named vector
# query <- KYCG_getDBs("KYCG.MM285.designGroup")[["TSS"]]
#
# result <- testEnrichmentSEA(query, "MM285.seqContextN")
# print(result[, c("dbname", "test", "estimate", "FDR", "nQ", "nD", "overlap")])
#
# # Plot by setting prepPlot=TRUE
# query <- KYCG_getDBs("KYCG.MM285.designGroup")[["TSS"]]
# db <- KYCG_getDBs("MM285.seqContextN", "distToTSS")
# res <- testEnrichmentSEA(query, db, prepPlot = TRUE)
# KYCG_plotSetEnrichment(res[[1]])
#
# # Alternatively, test enrichment with continuous query with dicrete database
# beta_values <- getBetas(sesameDataGet("MM285.1.SigDF"))
# res <- testEnrichmentSEA(beta_values, "MM285.chromHMM")
# res[, c("dbname", "test", "estimate", "FDR", "nQ", "nD", "overlap")]
#
# nolint end
