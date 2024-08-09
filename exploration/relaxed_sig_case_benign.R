library(stringr)
library(sesame)

din <- "/projects/p30791/methylation/differential_methylation"
dout <- "/projects/p30791/methylation/differential_methylation/KYCG"

dm_results <- list(
  OQ_vs_AN = read.csv(paste0(din, "/DML_results_refOQ_compAN.csv")),
  CUB_vs_OQ = read.csv(paste0(din, "/DML_results_refCUB_compOQ.csv"))
)

adj_p_thres <- 0.10
effect_thres <- 0.01

## Download public Databases
dbs <- KYCG_listDBGroups("EPIC")
dbs <- dbs[str_starts(dbs$Title, fixed("KYCG.EPIC.")), ]
db_list <- list()
for (db_id in dbs$Title) {
  db_list[[db_id]] <- KYCG_getDBs(db_id)
}

# tfbs_db <- KYCG_getDBs("KYCG.EPIC.TFBSconsensus.20211013")

## Get protein coding genes
txnranges <- sesameData_getTxnGRanges(sesameData_check_genome(NULL, platform = "EPIC"))
protein_coding_genes <- unique(
  data.frame(txnranges[txnranges$transcript_type == "protein_coding", ])[["gene_name"]]
)

for (setname in names(dm_results)) {
  print(paste0("#### ", setname, " ####"))
  # Identify the hyper and hypo probes
  df <- dm_results[[setname]]
  comp <- tail(strsplit(setname, "_")[[1]], 1) # comparison category

  hyper_results <- df[
    (p.adjust(df[[paste0("Pval_Sample.Region", comp)]], method = "BH") < adj_p_thres) &
      (df[[paste0("Est_Sample.Region", comp)]] >= effect_thres),
  ]
  hypo_results <- df[
    (p.adjust(df[[paste0("Pval_Sample.Region", comp)]], method = "BH") < adj_p_thres) &
      (df[[paste0("Est_Sample.Region", comp)]] <= (-1 * effect_thres)),
  ]

  probes <- list(
    hyper = hyper_results$Probe_ID,
    hypo = hypo_results$Probe_ID
  )

  print(paste0(
    length(probes[["hyper"]]), " probes hypermethylated, and ",
    length(probes$hypo), " probes hypomethylated."
  ))

  for (trend in names(probes)) {
    # Run enrichment analysis
    print(paste0("Testing ", setname, " ", trend, "methylation..."))
    query <- probes[[trend]]
    for (db_id in dbs$Title) {
      results <- testEnrichment(query, db_list[[db_id]], platform = "EPIC")

      # Apply q-value threshold
      results <- results[results$FDR < 0.05, ]
      # For ease of interpretation, we are only keeping enrichment and discarding depletion
      results <- results[results$estimate > 1.05, ]

      db_name <- strsplit(db_id, ".", fixed = TRUE)[[1]][3]

      if (dim(results)[1] > 0) {
        # Fix rows where p-value and FDR-adjusted p-values (q-value) suffer from underflow
        min_nonzero_value <- min(results$p.value[results$p.value > 0], na.rm = TRUE)
        results$p.value[results$p.value == 0] <- min_nonzero_value

        min_nonzero_value <- min(results$FDR[results$FDR > 0], na.rm = TRUE)
        results$FDR[results$FDR == 0] <- min_nonzero_value

        print(paste0(db_name, ": ", paste(results$dbname, collapse = ", ")))
        file_path <- paste0(
          dout, "/relaxed_testEnrichment_",
          trend, "_", setname, "_", db_name, ".csv"
        )
        write.csv(results, file_path, row.names = FALSE)
      } else {
        print(paste("Enrichment results for", db_name, "is empty!"))
      }
    }

    # Test enrichment for genes
    gene_dbs <- KYCG_buildGeneDBs(query, max_distance = 10000, platform = "EPIC")

    if (length(gene_dbs) > 0) {
      results <- testEnrichment(
        query,
        gene_dbs,
        platform = "EPIC"
      )
      results <- results[
        (results$FDR < 0.05) &
          (results$estimate > 1.05) &
          (results$gene_name %in% protein_coding_genes),
      ]

      if (dim(results)[1] > 0) {
        print(paste("Gene hits:", paste(results$gene_name, collapse = ", ")))
        file_path <- paste0(dout, "/relaxed_testEnrichment_", trend, "_", setname, "_genes.csv")
        write.csv(results, file_path, row.names = FALSE)
      } else {
        print("Gene results empty!")
      }
    } else {
      print("Gene DB empty!")
    }
  }
}
