library(stringr)
library(jsonlite)
library(sesame)
library(ggplot2)

probe_dir <- "/projects/p30791/methylation/differential_methylation"
dout <- probe_dir

## Read the query probe sets
probe_sets <- list(
  hyper_refHDB_compCUB = readLines(paste0(probe_dir, "/probe_set_hyper_refHDB_compCUB.txt")),
  hypo_refHDB_compCUB = readLines(paste0(probe_dir, "/probe_set_hypo_refHDB_compCUB.txt")),
  hyper_refAN_compTU = readLines(paste0(probe_dir, "/probe_set_hyper_refAN_compTU.txt")),
  hypo_refAN_compTU = readLines(paste0(probe_dir, "/probe_set_hypo_refAN_compTU.txt"))
)

## Download public Databases
db_id <- "KYCG.EPIC.chromHMM.20211020" # we're only using this one here
db <- KYCG_getDBs(db_id)

##########################################
#### Classify DM probes with chromHMM ####
##########################################

numbers <- as.numeric(sub("_.*", "", names(db)))
sorted_names <- names(db)[order(numbers)]

df <- data.frame(matrix(NA, nrow = length(sorted_names), ncol = 4))

rownames(df) <- sorted_names
colnames(df) <- c("HDBvsCUB.hyper", "HDBvsCUB.hypo", "ANvsTU.hyper", "ANvsTU.hypo")

ref_list <- c("HDB", "AN")
comp_list <- c("CUB", "TU")

for (chromname in sorted_names) {
  for (i in seq_along(ref_list)) {
    ref <- ref_list[i]
    comp <- comp_list[i]
    for (trend in c("hyper", "hypo")) {
      ct <- length(intersect(
        db[[chromname]],
        probe_sets[[paste0(trend, "_ref", ref, "_comp", comp)]]
      ))
      df[chromname, paste0(ref, "vs", comp, ".", trend)] <- ct
    }
  }
}

write.csv(
  df,
  file = paste0(dout, "/num_dm_probes_by_chromHMM.csv"),
  quote = FALSE,
  row.names = TRUE
)
