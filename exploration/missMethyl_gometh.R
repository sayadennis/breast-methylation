library(sesame)
library(minfi)
library(missMethyl)

din <- "/projects/p30791/methylation/differential_methylation"
dout <- din

## Read the query probe sets
probeset_names <- c(
  "hyper_ER-_refAN_compTU",
  "hyper_ER+_refAN_compTU",
  "hypo_ER-_refAN_compTU",
  "hypo_ER+_refAN_compTU",
  "hyper_refHDB_compCUB",
  "hypo_refHDB_compCUB",
  "hyper_refOQ_compAN",
  "hypo_refOQ_compAN",
  "hyper_refAN_compTU",
  "hypo_refAN_compTU"
)

probe_sets <- list()

for (setname in probeset_names) {
  file_path <- paste0(din, "/probe_set_", setname, ".txt")
  probe_sets[[setname]] <- readLines(file_path)
}

#############################################
#### First run gometh for all DML probes ####
#############################################

gometh_res <- list()

for (setname in probeset_names) {
  gometh_res[[setname]] <- list()
  test_go <- gometh(
    sig.cpg = probe_sets[[setname]],
    # all.cpg=rownames(top),  # TODO: output all CpGs used for DML during DML step!!
    collection = "GO",
    array.type = "EPIC",
    plot.bias = FALSE
  )
  test_kegg <- gometh(
    sig.cpg = probe_sets[[setname]],
    # all.cpg=rownames(top),  # TODO: output all CpGs used for DML during DML step!!
    collection = "KEGG",
    array.type = "EPIC",
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
  for (collection in c("GO", "KEGG")) {
    if (dim(gometh_res[[setname]][[collection]])[1] > 0) {
      write.csv(
        gometh_res[[setname]][[collection]],
        paste0(dout, "/gometh_", collection, "_", setname, ".csv"),
        row.names = TRUE
      )
    }
  }
}

######################################################
#### Next run gometh for probes driving TFBS hits ####
######################################################

## Download TFBS Database
db <- KYCG_getDBs("KYCG.EPIC.TFBSconsensus.20211013")

tfs <- c(
  "ZNF217", "AHR", "TLE3", "GREB1", "GRHL2",
  "ESR1", "PR", "FOSL2", "AR", "GATA3",
  "NR3C1", "HES2"
)

# Get the probe set names
left_name <- "hypo_refHDB_compCUB"
right_name <- "hypo_refAN_compTU"

# Get all probes in comparisons
left_probes <- probe_sets[[left_name]]
right_probes <- probe_sets[[right_name]]

# List of probes in the left- and right-comaprison DML sets
left_probes_tf_all <- c()
right_probes_tf_all <- c()

for (tf in tfs) {
  # Get probes that overlap with TFBS
  tf_probes <- db[[tf]]
  left_probes_tf <- intersect(left_probes, tf_probes)
  right_probes_tf <- intersect(right_probes, tf_probes)

  # Add to left and right probe sets that overlap with TFBS
  left_probes_tf_all <- union(left_probes_tf, left_probes_tf_all)
  right_probes_tf_all <- union(right_probes_tf, right_probes_tf_all)
}

# Get left- and right-only probes
tf_probes <- list(
  "TFBS_HDB_vs_CUB_only" = setdiff(left_probes_tf_all, right_probes_tf_all),
  "TFBS_AN_vs_TU_only" = setdiff(right_probes_tf_all, left_probes_tf_all),
  "TFBS_HDBvCUB_ANvTU_overlap" = intersect(left_probes_tf_all, right_probes_tf_all)
)

for (setname in names(tf_probes)) {
  for (collection in c("GO", "KEGG")) {
    # Run gometh enrichment analysis
    test <- gometh(
      sig.cpg = tf_probes[[setname]],
      collection = collection,
      array.type = "EPIC",
      plot.bias = FALSE
    )
    test <- test[test$P.DE < 0.05, ] # cut with nominal p-value in case no FDR hits
    if (dim(test)[1] > 0) {
      write.csv(
        test,
        paste0(dout, "/gometh_", collection, "_", setname, ".csv"),
        row.names = TRUE
      )
    }
  }
}
