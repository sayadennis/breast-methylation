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
betas = read.table("/projects/p30791/methylation/sesame_out/betas_filtered.csv", row.names=1, sep=",", header=TRUE)  # nrows=2000 for testing
meta = read.csv("/projects/p30791/methylation/data/meta.csv")

# TODO: look into why betas includes more samples names than meta??
betas = betas[,colnames(betas) %in% paste0("X", meta$IDAT)]
betas = as.matrix(betas)

## Exclude probes that are missing levels on variables of interest (sample region etc.)
cpg_ok = checkLevels(betas, meta$Sample.Region)
print(paste0(sum(cpg_ok), " probes have sufficient levels for sample region."))

betas = betas[cpg_ok,]

###########################################
#### Differential methylation analysis ####
###########################################

for (reference in c("Normal", "TU")) {
    ## Change metadata into factors
    meta$Sample.Region = relevel(factor(meta$Sample.Region), reference)
 
    ## Differential methylation analysis
    smry = DML(betas, ~Sample.Region + Age, meta=meta)  # Case.Control + BMI
    saveRDS(smry, paste0("/projects/p30791/methylation/sesame_out/model_summary_ref", reference, ".RDS"))
    test_result = summaryExtractTest(smry)
    write.csv(
        test_result, paste0("/projects/p30791/methylation/sesame_out/DML_results_ref", reference, ".csv"), 
        quote = FALSE, row.names = FALSE
    )
 
    ## DMR (Differentially methylated regions)
    for (contrast in dmContrasts(smry)) {
        merged = DMR(betas, smry, contrast, platform="EPIC") # merge probes to regions
        write.csv(
            merged, paste0("/projects/p30791/methylation/sesame_out/DMR_results_", contrast, "_ref", reference, ".csv"), 
            quote = FALSE, row.names = FALSE
        )
    }
}

df = data.frame(Age = meta$Age,
    BetaValue = betas[test_result$Probe_ID[nrow(test_result)],])

ggplot(df, aes(Age, BetaValue)) + geom_smooth(method="lm") + geom_point()
ggsave(paste0(plotdir, "/age_vs_betas.png"))

