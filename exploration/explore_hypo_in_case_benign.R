library(sesame)
library(dplyr)

din <- "/projects/p30791/methylation/differential_methylation/KYCG"

data <- read.csv(paste0(din, "/testEnrichment_Non_monotonic_valley_A_genes.csv"))

txnranges <- sesameData_getTxnGRanges(
  sesameData_check_genome(NULL, platform = "EPIC")
)
protein_coding_genes <- unique(
  data.frame(txnranges[txnranges$transcript_type == "protein_coding", ])[["gene_name"]]
)

hit_genes <- data$gene_name
protein_coding_hits <- intersect(hit_genes, protein_coding_genes)

writeLines(
  protein_coding_hits,
  paste0(din, "/Non_monotonic_valley_A_protein_coding_genes.txt")
)

cosmic <- read.csv("/projects/p30791/methylation/cosmic_census_genes.csv")

tsg <- cosmic %>%
  filter(grepl("TSG", Role.in.Cancer))

tsg_list <- c(unlist(strsplit(tsg$Synonyms, ",")), tsg$Gene.Symbol)
tsg_list <- unique(tsg_list)

hit_genes_tsg <- intersect(tsg_list, protein_coding_hits)

# > hit_genes_tsg
#  [1] "CARS1"    "ARHGAP26" "ARHGEF10" "BAP1"     "BUB1B"    "CIITA"
#  [7] "CTNNA1"   "DNM2"     "FANCC"    "FBXW7"    "FOXO3"    "HOXA11"
# [13] "HOXA9"    "NCOR2"    "NTHL1"    "PBRM1"    "RAD51B"   "SMARCE1"
# [19] "SPOP"     "STAT5B"   "TNFAIP3"  "TRAF7"    "TSC2"     "WNK2"
