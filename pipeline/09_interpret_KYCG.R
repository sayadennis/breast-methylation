library(stringr)
library(sesame)
library(ggplot2)
library(gridExtra)

din <- "/projects/p30791/methylation/differential_methylation"
dout <- din

## Read the query probe sets
probeset_names <- c(
  # Comparisons along TPX
  "hyper_refUN_compCUB",
  "hyper_refCUB_compOQ",
  "hyper_refOQ_compAN",
  "hyper_refAN_compTU",
  "hypo_refUN_compCUB",
  "hypo_refCUB_compOQ",
  "hypo_refOQ_compAN",
  "hypo_refAN_compTU",
  # Biomarker-related comparisons
  "hypo_ER-_refAN_compTU",
  "hyper_ER-_refAN_compTU",
  "hypo_ER+_refAN_compTU",
  "hyper_ER+_refAN_compTU"
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

db_id <- "KYCG.EPIC.TFBSconsensus.20211013"

db <- db_list[[db_id]]

###################################################
#### Collect CpG count and overlap information ####
###################################################

cols <- c(
  "left_trend", "right_trend", "left_comparison", "right_comparison", "TF",
  "left_only_size", "overlap_size", "right_only_size", "total", "left_only_pct",
  "overlap_pct", "right_only_pct"
)

cts <- data.frame(matrix(
  c(
    # hypo in CUB and TU
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "ZNF217", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "AHR", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "TLE3", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "GREB1", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "GRHL2", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "ESR1", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "PR", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "FOSL2", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "AR", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "GATA3", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "NR3C1", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "HES2", 0, 0, 0, 0, 0., 0., 0.,
    # hyper in AN and TU
    "hyper", "hyper", "refOQ_compAN", "refAN_compTU", "CDCA2", 0, 0, 0, 0, 0., 0., 0.,
    "hyper", "hyper", "refOQ_compAN", "refAN_compTU", "ATF7IP", 0, 0, 0, 0, 0., 0., 0.,
    "hyper", "hyper", "refOQ_compAN", "refAN_compTU", "ZFP57", 0, 0, 0, 0, 0., 0., 0.,
    # mixed
    "hyper", "hypo", "refOQ_compAN", "refOQ_compAN", "CNOT3", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hyper", "refUN_compCUB", "refAN_compTU", "TEAD1", 0, 0, 0, 0, 0., 0., 0.
  ),
  ncol = length(cols),
  byrow = TRUE
))
colnames(cts) <- cols

nrows <- dim(cts)[1]

for (i in 1:nrows) {
  # Get the probe set names
  left_name <- paste0(cts[i, "left_trend"], "_", cts[i, "left_comparison"])
  right_name <- paste0(cts[i, "right_trend"], "_", cts[i, "right_comparison"])
  # Get all probes in comparisons
  left_probes_all <- probe_sets[[left_name]]
  right_probes_all <- probe_sets[[right_name]]
  # Get probes that overlap with TFBS
  tf <- cts[i, "TF"]
  tf_probes <- db[[tf]]
  left_probes <- intersect(left_probes_all, tf_probes)
  right_probes <- intersect(right_probes_all, tf_probes)
  # count the set sizes
  left_only_size <- length(setdiff(left_probes, right_probes))
  right_only_size <- length(setdiff(right_probes, left_probes))
  overlap_size <- length(intersect(left_probes, right_probes))
  # Calculate percentages
  total <- length(unique(c(left_probes, right_probes)))
  left_only_pct <- 100 * left_only_size / total
  overlap_pct <- 100 * overlap_size / total
  right_only_pct <- 100 * right_only_size / total
  # Record to dataframe
  cts[i, "left_only_size"] <- left_only_size
  cts[i, "overlap_size"] <- overlap_size
  cts[i, "right_only_size"] <- right_only_size
  cts[i, "total"] <- total
  cts[i, "left_only_pct"] <- left_only_pct
  cts[i, "overlap_pct"] <- overlap_pct
  cts[i, "right_only_pct"] <- right_only_pct
}

write.csv(cts, file = paste0(dout, "/TFBS_hits_probe_overlaps.csv"), row.names = FALSE)

####################################################
#### Collect gene count and overlap information ####
####################################################

cts <- data.frame(matrix(
  c(
    # hypo in CUB and TU
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "ZNF217", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "AHR", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "TLE3", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "GREB1", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "GRHL2", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "ESR1", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "PR", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "FOSL2", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "AR", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "GATA3", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "NR3C1", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hypo", "refUN_compCUB", "refAN_compTU", "HES2", 0, 0, 0, 0, 0., 0., 0.,
    # hyper in AN and TU
    "hyper", "hyper", "refOQ_compAN", "refAN_compTU", "CDCA2", 0, 0, 0, 0, 0., 0., 0.,
    "hyper", "hyper", "refOQ_compAN", "refAN_compTU", "ATF7IP", 0, 0, 0, 0, 0., 0., 0.,
    "hyper", "hyper", "refOQ_compAN", "refAN_compTU", "ZFP57", 0, 0, 0, 0, 0., 0., 0.,
    # mixed
    "hyper", "hypo", "refOQ_compAN", "refOQ_compAN", "CNOT3", 0, 0, 0, 0, 0., 0., 0.,
    "hypo", "hyper", "refUN_compCUB", "refAN_compTU", "TEAD1", 0, 0, 0, 0, 0., 0., 0.
  ),
  ncol = length(cols),
  byrow = TRUE
))
colnames(cts) <- cols

nrows <- dim(cts)[1]

txnranges <- sesameData_getTxnGRanges(sesameData_check_genome(NULL, platform = "EPIC"))
protein_coding_genes <- unique(
  data.frame(txnranges[txnranges$transcript_type == "protein_coding", ])[["gene_name"]]
)

for (i in 1:nrows) {
  # Get the probe set names
  left_name <- paste0(cts[i, "left_trend"], "_", cts[i, "left_comparison"])
  right_name <- paste0(cts[i, "right_trend"], "_", cts[i, "right_comparison"])
  # Get all probes in comparisons
  left_probes_all <- probe_sets[[left_name]]
  right_probes_all <- probe_sets[[right_name]]
  # Get probes that overlap with TFBS
  tf <- cts[i, "TF"]
  tf_probes <- db[[tf]]
  left_probes <- intersect(left_probes_all, tf_probes)
  right_probes <- intersect(right_probes_all, tf_probes)
  # Get gene lists
  left_res <- testEnrichmentGene(left_probes, platform = "EPIC")
  right_res <- testEnrichmentGene(right_probes, platform = "EPIC")
  left_sig_genes <- left_res[left_res$FDR < 0.05, "gene_name"]
  right_sig_genes <- right_res[right_res$FDR < 0.05, "gene_name"]
  left_sig_genes <- intersect(left_sig_genes, protein_coding_genes)
  right_sig_genes <- intersect(right_sig_genes, protein_coding_genes)
  # get unique and overlapping genes
  left_only_genes <- setdiff(left_sig_genes, right_sig_genes)
  right_only_genes <- setdiff(right_sig_genes, left_sig_genes)
  overlap_genes <- intersect(left_sig_genes, right_sig_genes)
  # count the set sizes
  left_only_size <- length(left_only_genes)
  right_only_size <- length(right_only_genes)
  overlap_size <- length(overlap_genes)
  # Calculate percentages
  total <- length(unique(c(left_sig_genes, right_sig_genes)))
  left_only_pct <- 100 * left_only_size / total
  overlap_pct <- 100 * overlap_size / total
  right_only_pct <- 100 * right_only_size / total
  # Record to dataframe
  cts[i, "left_only_size"] <- left_only_size
  cts[i, "overlap_size"] <- overlap_size
  cts[i, "right_only_size"] <- right_only_size
  cts[i, "total"] <- total
  cts[i, "left_only_pct"] <- left_only_pct
  cts[i, "overlap_pct"] <- overlap_pct
  cts[i, "right_only_pct"] <- right_only_pct
  # write genes to TXT
  writeLines(left_only_genes, paste0(dout, "/genes_", tf, "_hits_", left_name, "_only.txt"))
  writeLines(right_only_genes, paste0(dout, "/genes_", tf, "_hits_", right_name, "_only.txt"))
  writeLines(
    overlap_genes,
    paste0(dout, "/genes_overlap_", tf, "_hits_", left_name, "_AND_", right_name, "_only.txt")
  )
}

## Get genes targeted by PRC2 components of hypermethylation TFBS hits
for (tf in c("JARID2", "TCF21", "EZH2", "SUZ12")) {
  probes_all <- probe_sets[["hyper_refAN_compTU"]]
  tf_probes <- db[[tf]]
  probes <- intersect(probes_all, tf_probes)
  res <- testEnrichmentGene(probes, platform = "EPIC")
  sig_genes <- res[res$FDR < 0.05, "gene_name"]
  sig_genes <- intersect(sig_genes, protein_coding_genes)
  writeLines(sig_genes, paste0(dout, "/genes_", tf, "_hits_hyper_refAN_compTU.txt"))
}

write.csv(cts, file = paste0(dout, "/TFBS_hits_gene_overlaps.csv"), row.names = FALSE)
