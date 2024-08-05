library(rjson)
library(sesame)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)

plotdir <- "/projects/p30791/methylation/plots"
raw_dir <- "/projects/p30791/methylation/raw_data"
sesame_dir <- "/projects/p30791/methylation/sesame_data"
dout <- "/projects/p30791/methylation/differential_methylation"

############################
#### Load and prep data ####
############################

config <- fromJSON(file = "~/breast-methylation/pipeline/config.json")
p_thres <- config$p_thres
effect_thres <- config$effect_thres

## Load data
betas <- read.table(
  paste0(sesame_dir, "/betas_processed.csv"),
  row.names = 1, sep = ",", header = TRUE
) # nrows=2000 for testing
meta <- read.csv(paste0(raw_dir, "/meta.csv"))

betas <- as.matrix(betas)
meta <- meta[paste0("X", meta$IDAT) %in% colnames(betas), ]

# Mean imputation for BMI's missing values (control during DM analysis)
meta$BMI <- ifelse(is.na(meta$BMI), mean(meta$BMI, na.rm = TRUE), meta$BMI)

##########################
#### Define functions ####
##########################

get_probe_sets <- function(test_result, p_thres, effect_thres, pval_colname, slope_colname) {
  # Define boolean vector for significance
  pvals <- test_result[[pval_colname]]
  pvals[is.na(pvals)] <- 0.99 # fill NA with non-significant number
  pvals[pvals == 0] <- min(pvals[pvals != 0]) # fix underflow
  significant_bool <- p.adjust(pvals, method = "BH") < p_thres

  # Get hyper and hypo methylated probe IDs
  significant_hyper <- test_result[
    significant_bool & (test_result[[slope_colname]] >= effect_thres), "Probe_ID"
  ]
  significant_hypo <- test_result[
    significant_bool & (test_result[[slope_colname]] <= (-1) * effect_thres), "Probe_ID"
  ]

  # Convert to characterstring vectors
  significant_hyper <- as.character(significant_hyper[[1]])
  significant_hypo <- as.character(significant_hypo[[1]])

  # Return as named list
  list(
    hyper = significant_hyper,
    hypo = significant_hypo
  )
}

##############################################################
#### Differential methylation analysis by tissue category ####
##############################################################

refs <- c("CFN", "CUB", "OQ", "AN", "CFN", "CFN", "CUB", "OQ", "AN")
comps <- c("TU", "TU", "TU", "TU", "AN", "CUB", "OQ", "AN", "TU")

for (i in seq_along(refs)) {
  ref <- refs[i]
  comp <- comps[i]

  ## Subset the betas and metadata
  meta_sub <- meta[meta$Sample.Region %in% c(ref, comp), ]
  meta_sub$Sample.Region <- relevel(factor(meta_sub$Sample.Region), ref)
  betas_sub <- betas[, paste0("X", meta_sub$IDAT)]

  ## Exclude probes that are missing levels on sample region etc.
  cpg_ok <- checkLevels(betas_sub, meta_sub$Sample.Region)
  print(paste0(sum(cpg_ok), " probes have sufficient levels for sample region."))
  betas_sub <- betas_sub[cpg_ok, ]

  ## Differential methylation analysis
  smry <- DML(betas_sub, as.formula("~ Sample.Region + Age + BMI"), meta = meta_sub)
  saveRDS(
    smry,
    paste0(
      dout, "/model_summary_ref", ref, "_comp", comp, ".RDS"
    )
  )
  test_result <- summaryExtractTest(smry)
  # Save all results
  write.csv(
    test_result,
    paste0(
      dout, "/DML_results_ref", ref, "_comp", comp, ".csv"
    ),
    quote = FALSE, row.names = FALSE
  )

  # Save significant probes
  pval_colname <- paste0("Pval_Sample.Region", comp)
  slope_colname <- paste0("Est_Sample.Region", comp)

  probe_sets <- get_probe_sets(
    test_result = test_result,
    p_thres = p_thres,
    effect_thres = effect_thres,
    pval_colname = pval_colname,
    slope_colname = slope_colname
  )

  writeLines(
    probe_sets$hyper,
    paste0(dout, "/probe_set_hyper_ref", ref, "_comp", comp, ".txt")
  )
  writeLines(
    probe_sets$hypo,
    paste0(dout, "/probe_set_hypo_ref", ref, "_comp", comp, ".txt")
  )

  ## DMR (Differentially methylated regions)
  for (contrast in dmContrasts(smry)) {
    # merge probes to regions
    merged <- DMR(betas_sub, smry, contrast, platform = "EPIC")
    write.csv(
      merged,
      paste0(
        dout, "/DMR_results_", contrast, "_ref", ref, ".csv"
      ),
      quote = FALSE, row.names = FALSE
    )
  }
}

# nolint start
# df <- data.frame(
#   Age = meta$Age,
#   BetaValue = betas[test_result$Probe_ID[nrow(test_result)], ]
# )
#
# ggplot(df, aes(Age, BetaValue)) +
#   geom_smooth(method = "lm") +
#   geom_point()
# ggsave(paste0(plotdir, "/age_vs_betas.png"))
# nolint end

###################################################################
#### Differential methylation along TPX separated by ER status ####
###################################################################

refs <- c("CUB", "OQ", "AN")
comps <- c("OQ", "AN", "TU")

for (i in seq_along(refs)) {
  ref <- refs[i]
  comp <- comps[i]

  # Subset the betas and metadata
  meta_sub <- meta[meta$Sample.Region %in% c(ref, comp), ]
  betas_sub <- betas[, paste0("X", meta_sub$IDAT)]

  for (er_status in unique(meta_sub$ER)) {
    # Subset metadata and betas by ER status
    meta_sub_er <- meta_sub[meta_sub$ER == er_status, ]
    meta_sub_er$Sample.Region <- relevel(factor(meta_sub_er$Sample.Region), ref)
    betas_sub_er <- betas_sub[, paste0("X", meta_sub_er$IDAT)]

    # Exclude probes that are missing levels on sample region etc.
    cpg_ok <- checkLevels(betas_sub_er, meta_sub_er$Sample.Region)
    print(paste0(sum(cpg_ok), " probes have sufficient levels for sample region."))
    betas_sub_cpgfilter <- betas_sub_er[cpg_ok, ]

    # Run DML
    smry <- DML(betas_sub_cpgfilter, as.formula("~ Sample.Region + Age + BMI"), meta = meta_sub_er)
    test_result <- summaryExtractTest(smry)

    # Save results
    write.csv(
      test_result,
      paste0(
        dout, "/DML_results_ER", er_status, "_only_", ref, "_vs_", comp, ".csv"
      ),
      quote = FALSE, row.names = FALSE
    )

    # Save significant probes
    pval_colname <- paste0("Pval_Sample.Region", comp)
    slope_colname <- paste0("Est_Sample.Region", comp)

    probe_sets <- get_probe_sets(
      test_result = test_result,
      p_thres = p_thres,
      effect_thres = effect_thres,
      pval_colname = pval_colname,
      slope_colname = slope_colname
    )

    writeLines(
      probe_sets$hyper,
      paste0(dout, "/probe_set_hyper_ER", er_status, "_ref", ref, "_comp", comp, ".txt")
    )
    writeLines(
      probe_sets$hypo,
      paste0(dout, "/probe_set_hypo_ER", er_status, "_ref", ref, "_comp", comp, ".txt")
    )
  }
}

#########################################################################################
#### Differential methylation between same tissue types separated by tumor biomarker ####
#########################################################################################

comps <- c("CUB", "OQ", "AN", "TU")

meta_cases <- meta[meta$Case.Control == "case", ]
betas_cases <- betas[, paste0("X", meta_cases$IDAT)]

# Create useful metadata columns
meta_cases$HER2 <- relevel(
  factor(
    ifelse(
      meta_cases$HER2 %in% c(2, 3),
      "Positive", "Negative"
    )
  ), "Negative"
)
meta_cases$ER <- relevel(
  factor(
    ifelse(
      meta_cases$ER == "+",
      "Positive", "Negative"
    )
  ), "Negative"
)

for (tissue_category in comps) {
  print(paste0("Running analysis for ", tissue_category, "..."))
  # Subset the betas and metadata
  meta_sub <- meta_cases[meta_cases$Sample.Region == tissue_category, ]
  betas_sub <- betas_cases[, paste0("X", meta_sub$IDAT)]

  for (contrast_feature in c("ER", "HER2")) {
    print(paste0("Running ", contrast_feature, "..."))
    # Exclude probes that are missing levels on sample region etc.
    cpg_ok <- checkLevels(betas_sub, meta_sub[[contrast_feature]])
    print(paste0(sum(cpg_ok), " probes have sufficient levels for sample region."))
    betas_sub_cpgfilter <- betas_sub[cpg_ok, ]
    # Create formula string with contrast feature to plug into DML
    formula_string <- paste("~", contrast_feature, "+ Age + BMI")
    formula_expression <- as.formula(formula_string)
    # Run DML
    smry <- DML(betas_sub_cpgfilter, formula_expression, meta = meta_sub)
    test_result <- summaryExtractTest(smry)
    # Save results
    write.csv(
      test_result,
      paste0(
        dout, "/DML_results_of_", tissue_category, "_by_", contrast_feature, ".csv"
      ),
      quote = FALSE, row.names = FALSE
    )
    # Save significant probes
    pval_colname <- paste0("Pval_", contrast_feature, "Positive")
    slope_colname <- paste0("Est_", contrast_feature, "Positive")

    probe_sets <- get_probe_sets(
      test_result = test_result,
      p_thres = p_thres,
      effect_thres = effect_thres,
      pval_colname = pval_colname,
      slope_colname = slope_colname
    )

    # Write probe sets
    writeLines(
      probe_sets$hyper,
      paste0(
        dout,
        "/probe_set_hyper_",
        contrast_feature,
        "_neg_vs_pos_in_",
        tissue_category,
        ".txt"
      )
    )
    writeLines(
      probe_sets$hypo,
      paste0(
        dout,
        "/probe_set_hyper_",
        contrast_feature,
        "_neg_vs_pos_in_",
        tissue_category,
        ".txt"
      )
    )
  }
}
