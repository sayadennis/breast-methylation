library(stringr)
library(sesame)

din <- "/projects/p30791/methylation/differential_methylation"

dm_results <- list(
  OQ_vs_AN = read.csv(paste0(din, "/DML_results_refOQ_compAN.csv")),
  CUB_vs_OQ = read.csv(paste0(din, "/DML_results_refCUB_compOQ.csv"))
)

nominal_p_thres <- 0.05
effect_thres <- 0.05

## Download public Databases
tfbs_db <- KYCG_getDBs("KYCG.EPIC.TFBSconsensus.20211013")

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
    (df[[paste0("Pval_Sample.Region", comp)]] < nominal_p_thres) &
      (df[[paste0("Est_Sample.Region", comp)]] >= effect_thres),
  ]
  hypo_results <- df[
    (df[[paste0("Pval_Sample.Region", comp)]] < nominal_p_thres) &
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
    results <- testEnrichment(query, tfbs_db, platform = "EPIC")

    # Apply q-value threshold
    results <- results[results$FDR < 0.05, ]
    # For ease of interpretation, we are only keeping enrichment and discarding depletion
    results <- results[results$estimate > 1.05, ]

    if (dim(results)[1] == 0) {
      # only save results if any are significant
      print("TFBS results empty!")
      next
    }

    # Fix rows where p-value and FDR-adjusted p-values (q-value) suffer from underflow
    min_nonzero_value <- min(results$p.value[results$p.value > 0], na.rm = TRUE)
    results$p.value[results$p.value == 0] <- min_nonzero_value

    min_nonzero_value <- min(results$FDR[results$FDR > 0], na.rm = TRUE)
    results$FDR[results$FDR == 0] <- min_nonzero_value

    print(paste("TFBS hits:", paste(results$dbname, collapse = ", ")))

    # Test enrichment for genes
    gene_dbs <- KYCG_buildGeneDBs(query, max_distance = 10000, platform = "EPIC")

    if (length(gene_dbs) == 0) {
      print("Gene DB empty!")
      next
    }

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

    if (dim(results)[1] == 0) {
      print("Gene results empty!")
      next
    }

    print(paste("Gene hits:", paste(results$gene_name, collapse = ", ")))
  }
}
