library(dplyr)
library(tidyr)
library(sesame)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Please provide 3 arguments: 1) data (?) directory, 2) metadata file path, and 3) output to write predicted meta information.", call.=FALSE)
} else {
  ddata = args[1]  # e.g. "/projects/p30791/methylation/sesame_out"
  meta_fn = args[2]  # e.g. "/projects/p30791/methylation/data/meta.csv"
  dout = args[3]  # e.g. "/projects/p30791/methylation/sesame_out/data_summary"
}

## download clock file from http://zwdzwd.github.io/InfiniumAnnotation
#model <- readRDS("~/Downloads/Clock_Horvath353.rds")
################################################################

meta = read.csv(meta_fn)
sdfs = readRDS(paste0(ddata, "/sdf_processed.RDS"))
betas = read.csv(paste0(ddata, "/betas_processed.csv"), row.names=1)

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

write.table(meta, paste0(dout,"/predicted_meta.csv"), row.names = FALSE, quote = FALSE, sep=",")

female_ratio <- sum(meta$InferredSex == 'FEMALE') / nrow(meta)
cat(sprintf("Value counts of meta sex: %.1f%% FEMALE\n", 100 * female_ratio))

cts <- table(meta$Race, meta$InferredEthnicity)
print(cts)

