library(dplyr)
library(tidyr)
library(sesame)

ddata = "/projects/p30791/methylation/sesame_out"
dmeta = "/projects/p30791/methylation/data"
dout = "/projects/p30791/methylation/plots"

## download clock file from http://zwdzwd.github.io/InfiniumAnnotation
#model <- readRDS("~/Downloads/Clock_Horvath353.rds")
################################################################

meta = read.csv(paste0(dmeta, "/meta.csv"))
sdfs = readRDS(paste0(ddata, "/sdf_prepped.RDS"))
betas = read.csv(paste0(ddata, "/betas.csv"), row.names=1)

for (idat_id in names(sdfs)) {
    sdf = sdfs[idat_id][[1]]
    meta[meta$IDAT==idat_id,"InferredSex"] = inferSex(sdf)
    meta[meta$IDAT==idat_id,"InferredSexKaryotypes"] = inferSexKaryotypes(sdf)
    meta[meta$IDAT==idat_id,"InferredEthnicity"] = inferEthnicity(sdf)
    if (idat_id %in% colnames(betas)) {
        betas_sample <- betas[idat_id]
        names(betas_sample) <- rownames(betas)
        meta[meta$IDAT==idat_id,"EstLeukFrac"] = estimateLeukocyte(betas_sample)
    } else {
        print(paste0("IDAT ", idat_id, " was not in the betas CSV columns."))
    }
    #meta[meta$IDAT==idat_id,"PredictedAge"] = predictAge(betas_sample, model)  # cannot find model in specified link
}

write.table(meta, "/home/srd6051/predicted_meta.csv", row.names = FALSE, quote = FALSE, sep=",")

female_ratio <- sum(meta$InferredSex == 'FEMALE') / nrow(meta)
cat(sprintf("Value counts of meta sex: %.1f%% FEMALE\n", 100 * female_ratio))

cts <- table(meta$Race, meta$InferredEthnicity)
print(cts)

