library(sesame)
library(minfi)
library(missMethyl)

din <- "/projects/p30791/methylation/sesame_out/differential_methylation"
dout <- din

## Read the query probe sets
probeset_names <- c(
  "hypo_in_TU_ER",
  "hyper_in_TU_ER",
  "CUB_down_from_CFN",
  "CUB_up_from_CFN",
  "AN_up_from_OQ",
  "AN_down_from_OQ",
  "TU_down_from_AN",
  "TU_up_from_AN",
  "Monotonic_increase_A",
  "Monotonic_increase_B",
  "Monotonic_increase_C",
  "Monotonic_decrease_A",
  "Monotonic_decrease_B",
  "Monotonic_decrease_C"
)

probe_sets <- list()

for (setname in probeset_names) {
  file_path <- paste0(din, "/probe_set_", setname, ".txt")
  probe_sets[[setname]] <- readLines(file_path)
}

gometh_res <- list()

for (setname in probeset_names) {
  gometh_res[[setname]] <- list()
  test_go <- gometh(
    sig.cpg = probe_sets[[setname]],
    # all.cpg=rownames(top),  # TODO: output all CpGs used for DML during DML step!!
    collection = c("GO", "KEGG"),
    plot.bias = FALSE
  )
  test_kegg <- gometh(
    sig.cpg = probe_sets[[setname]],
    # all.cpg=rownames(top),  # TODO: output all CpGs used for DML during DML step!!
    collection = "KEGG",
    plot.bias = FALSE
  )
  test_go <- test_go[test_go$FDR < 0.05, ]
  test_kegg <- test_kegg[test_kegg$FDR < 0.05, ]
  gometh_res[[setname]][["GO"]] <- test_go
  gometh_res[[setname]][["KEGG"]] <- test_kegg
  print(paste("####", setname, "####"))
  print(paste("Number of CpGs in set:", length(probe_sets[[setname]])[1]))
  print(paste("TOTAL HITS:", dim(test_go)[1] + dim(test_kegg)[1]))
  print(paste("  GO:", dim(test_go)[1]))
  print(paste("  KEGG:", dim(test_kegg)[1]))
  print("")
  # Write results to disk
  for (query in c("GO", "KEGG")) {
    if (dim(gometh_res[[setname]][[query]])[1] > 0) {
      write.csv(
        gometh_res[[setname]][[query]],
        paste0(dout, "/gometh_", query, "_", setname, ".csv"),
        row.names = TRUE
      )
    }
  }
}
