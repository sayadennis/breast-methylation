library(sesame)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)

plotdir = "/projects/p30791/methylation/plots"

############################
#### Load and prep data ####
############################

## Load data
betas = read.table("/projects/p30791/methylation/sesame_out/betas.csv", row.names=1, sep=",", header=TRUE)  # nrows=2000 for testing
meta = read.csv("/projects/p30791/methylation/data/meta.csv")

# TODO: look into why betas includes more samples names than meta??
betas = betas[,colnames(betas) %in% paste0("X", meta$IDAT)]
betas = as.matrix(betas)

## Exclude probes that are missing levels on variables of interest (sample region etc.)
cpg_ok = (
    checkLevels(betas, meta$Sample.Region) &
    checkLevels(betas, meta$Case.Control)
)
print(paste0(sum(cpg_ok), " probes have sufficient levels for sample region."))

betas = betas[cpg_ok,]

## Change metadata into factors
meta$Sample.Region = relevel(factor(meta$Sample.Region), "Normal")
meta$Case.Control = relevel(factor(meta$Case.Control), "normal")

###########################################
#### Differential methylation analysis ####
###########################################

## Differential methylation analysis
smry = DML(betas, ~Sample.Region + Age, meta=meta)  # Case.Control + BMI

saveRDS(smry, "/projects/p30791/methylation/sesame_out/model_summary.RDS")

test_result = summaryExtractTest(smry)

test = test_result %>% 
    dplyr::filter(FPval_Sample.Region < 0.05, Eff_Sample.Region > 0.1) %>% 
    select(FPval_Sample.Region, Eff_Sample.Region)

test_result %>% dplyr::select(Probe_ID, Est_Age, Pval_Age) %>% tail

df = data.frame(Age = meta$Age,
    BetaValue = betas[test_result$Probe_ID[nrow(test_result)],])

ggplot(df, aes(Age, BetaValue)) + geom_smooth(method="lm") + geom_point()
ggsave(paste0(plotdir, "/age_vs_betas.png"))

## DMR (Differentially methylated regions)
dmContrasts(smry)  # pick a contrast from output

merged = DMR(betas, smry, "Sample.RegionTU", platform="EPIC") # merge probes to regions

sig = merged %>% dplyr::filter(Seg_Pval_adj < 0.01)  # statistically significant probes and regions

