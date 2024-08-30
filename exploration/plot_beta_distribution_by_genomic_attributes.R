library(stringr)
library(sesame)
library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)

plotdir <- "/projects/p30791/methylation/plots"
raw_dir <- "/projects/p30791/methylation/raw_data"
sesame_dir <- "/projects/p30791/methylation/sesame_data"

## Load betas and meta data
betas <- read.table(
  paste0(sesame_dir, "/betas_processed.csv"),
  row.names = 1, sep = ",", header = TRUE
) # nrows=2000 for testing
meta <- read.csv(paste0(raw_dir, "/meta.csv"))

meta <- meta[paste0("X", meta$IDAT) %in% colnames(betas), ]

# Mean imputation for BMI's missing values (control during DM analysis)
meta$BMI <- ifelse(is.na(meta$BMI), mean(meta$BMI, na.rm = TRUE), meta$BMI)

## Download public Databases
dbs <- KYCG_listDBGroups("EPIC")
dbs <- dbs[str_starts(dbs$Title, fixed("KYCG.EPIC.")), ]
db_list <- list()
for (db_id in dbs$Title) {
  db_list[[db_id]] <- KYCG_getDBs(db_id)
}

## Process data into long form for violin plots
# convert to long form
betas_long <- betas %>%
  rownames_to_column(var = "Probe_ID") %>%
  pivot_longer(cols = -Probe_ID, names_to = "Sample", values_to = "Beta")

# add tissue category information to long form
betas_long <- betas_long %>%
  mutate(Sample_ID = gsub("^X", "", Sample)) %>%
  left_join(
    meta %>%
      select(IDAT, Sample.Region),
    by = c("Sample_ID" = "IDAT")
  ) %>%
  select(-Sample_ID)

# # add CGI metadata
# epic <- read.csv(
#   "/projects/p30791/methylation/illumina_array_meta/EPIC-8v2-0_A1.csv",
#   skip=7, nrows=10000)
# epic <- epic[c("Name", "Relation_to_UCSC_CpG_Island")]
#
# epic_long <- epic %>%
#   separate_rows(Relation_to_UCSC_CpG_Island, sep = ";")
#
# merged_df <- betas_long %>%
#   left_join(epic_long, by = c("Probe_ID" = "Name"), relationship = "many-to-many")
#
# merged_df <- drop_na(merged_df)
# merged_df <- filter(merged_df, Relation_to_UCSC_CpG_Island != "")
#
# mean_betas <- merged_df %>%
#   group_by(Sample, Sample.Region, Relation_to_UCSC_CpG_Island) %>%
#   summarise(Beta = mean(Beta, na.rm = TRUE), .groups = 'drop')

##################
#### Plot CGI ####
##################

# Use SeSAMe CGI database instead
db <- db_list[["KYCG.EPIC.CGI.20210713"]]
all_probes <- unique(unlist(db))
na_matrix <- matrix(NA, nrow = length(all_probes), ncol = 1)
sesame_cgi <- as.data.frame(na_matrix)
rownames(sesame_cgi) <- all_probes
colnames(sesame_cgi) <- c("CGI_location")

for (locname in c("Island", "N_Shelf", "N_Shore", "S_Shelf", "S_Shore")) {
  sesame_cgi[db[[locname]], "CGI_location"] <- locname
}

sesame_cgi <- rownames_to_column(sesame_cgi, var = "Probe_ID")

merged_df_sesame <- betas_long %>%
  left_join(sesame_cgi, by = c("Probe_ID" = "Probe_ID"))

merged_df_sesame <- drop_na(merged_df_sesame)

mean_betas <- merged_df_sesame %>%
  group_by(Sample, Sample.Region, CGI_location) %>%
  summarise(Beta = mean(Beta, na.rm = TRUE), .groups = "drop")

## Plot

mean_betas$Sample.Region <- factor(
  mean_betas$Sample.Region,
  levels = c("CFN", "CUB", "OQ", "AN", "TU")
)
mean_betas$CGI_location <- factor(
  mean_betas$CGI_location,
  levels = c("N_Shelf", "N_Shore", "Island", "S_Shore", "S_Shelf")
)
colors <- c(
  "CFN" = "gray66",
  "CUB" = "steelblue2",
  "OQ" = "chartreuse3",
  "AN" = "orange1",
  "TU" = "red"
)

p <- ggplot(mean_betas, aes(x = CGI_location, y = Beta, fill = Sample.Region)) +
  # geom_boxplot(outlier.size = 0.15, size = 0.15) +
  geom_violin(linewidth = 0.4) +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "Beta values", fill = "Tissue Category")

ggsave(paste0(plotdir, "/cgi_violinplot.png"), plot = p, width = 8, height = 2, units = "in")

write.csv(mean_betas, paste0(plotdir, "/mean_long_betas_cgi.csv"))

#######################
#### Plot chromHMM ####
#######################

# Use SeSAMe chromHMM database instead
db <- db_list[["KYCG.EPIC.chromHMM.20211020"]]
all_probes <- unique(unlist(db))
na_matrix <- matrix(NA, nrow = length(all_probes), ncol = 1)
sesame_chromhmm <- as.data.frame(na_matrix)
rownames(sesame_chromhmm) <- all_probes
colnames(sesame_chromhmm) <- c("chromHMM")

ordered_names <- names(db)[order(as.numeric(sapply(strsplit(names(db), "_"), `[`, 1)))]

for (locname in names(db)) {
  sesame_chromhmm[db[[locname]], "chromHMM"] <- locname
}

sesame_chromhmm <- rownames_to_column(sesame_chromhmm, var = "Probe_ID")

merged_df_sesame <- betas_long %>%
  left_join(sesame_chromhmm, by = c("Probe_ID" = "Probe_ID"))

merged_df_sesame <- drop_na(merged_df_sesame)

mean_betas <- merged_df_sesame %>%
  group_by(Sample, Sample.Region, chromHMM) %>%
  summarise(Beta = mean(Beta, na.rm = TRUE), .groups = "drop")

## Plot

mean_betas$Sample.Region <- factor(
  mean_betas$Sample.Region,
  levels = c("CFN", "CUB", "OQ", "AN", "TU")
)
mean_betas$chromHMM <- factor(
  mean_betas$chromHMM,
  levels = ordered_names
)
colors <- c(
  "CFN" = "gray66",
  "CUB" = "steelblue2",
  "OQ" = "chartreuse3",
  "AN" = "orange1",
  "TU" = "red"
)

p <- ggplot(mean_betas, aes(x = chromHMM, y = Beta, fill = Sample.Region)) +
  geom_boxplot(outlier.size = 0.05, size = 0.15) +
  # geom_violin(linewidth=0.4) +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "Beta values", fill = "Tissue Category") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(paste0(plotdir, "/chromHMM_boxplot.png"), plot = p, width = 16, height = 3, units = "in")
