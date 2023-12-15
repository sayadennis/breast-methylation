library(stringr)
library(sesame)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Please provide 3 arguments: 1) input directory, 2) output directory to write plots to, and 3) metadata file path.", call.=FALSE)
} else {
  din = args[1]  # e.g. "/projects/p30791/methylation/sesame_out"
  dout = args[2]  # e.g. "/projects/p30791/methylation/plots"
  meta_fn = args[3]  # e.g. "/projects/p30791/methylation/data/meta.csv"
}

## Load data
sdfs = readRDS(paste0(din, "/sdf_prepped.RDS"))
qcs = readRDS(paste0(din, "/qcs.RDS"))
qc_df = read.csv(paste0(din, "/qc_metrics.csv"))
meta = read.csv(meta_fn)

sesameQC_plotBar_custom <- function(qcs, key) {
    if (is(qcs, "sesameQC")) { qcs <- list(qcs); }
    df <- do.call(rbind, lapply(qcs, function(x) as.data.frame(x@stat)))

    ## set display names
    g <- qcs[[1]]@group
    display_nms <- do.call(c, lapply(names(g), function(gn) { setNames(
        sprintf("%s | %s", gn, str_trim(g[[gn]])), names(g[[gn]]))}))

    if (ncol(df) == 0) { stop("There is no QC metrics to plot") }

    df$sample_name <- names(qcs) # infer sample name
    
    plt <- NULL
    x <- key
    p <- ggplot(df) +
        geom_bar(aes_string("sample_name", x), stat="identity") +
        ylab("") + xlab("") + ggtitle(display_nms[x]) +
        theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
    ## customization for the most important
    if (x == "frac_dt") {
        p <- p + scale_y_continuous(labels=scales::percent)
    } else if (x == "median_beta_cg" || x == "median_beta_ch") {
        p <- p + ylim(c(0,1))
    }
    plt
}

## Plot QC bar graphs
for (groupname in c("TU", "AN", "OQ", "CUB", "Normal")) {
    for (key in c("frac_dt", "mean_intensity", "median_beta_cg", "median_beta_ch", "RGratio", "RGdistort")) {
        if (groupname=="Normal") {
            ratio = 2.5
        } else {
            ratio = 1.0
        }
        qc_subset = qcs[meta[meta["Sample.Region"]==groupname,]$IDAT]
        png(filename = paste0(dout, "/qc_plotbar_", key, "_", groupname, ".png"), width=ratio * 7.5, height=1.0 * 7.5, type="cairo")
        sesameQC_plotBar_custom(qc_subset, key=key)
        dev.off()
        #ggsave(filename=paste0(dout, "/qc_plotbar_", key, "_", groupname, ".png"), width=ratio * 7.5, height=1.0 * 7.5)
    }
}

# Note samples with high distortion
distorted_sampleids = qc_df[abs(qc_df$RGdistort-1)>0.5,]$X
print(paste0("Samples with >0.5 distortion: ", paste0(distorted_sampleids, collapse=", ")))
print("Meta data of these samples:")
print(meta[meta$IDAT %in% distorted_sampleids,])
writeLines(distorted_sampleids, paste0(din, "/exclude_IDATs.txt"))

## Plot QC scatter plots
plot_sampleids = c(distorted_sampleids, meta$IDAT[1])
for (idat_id in plot_sampleids) {
    sdf = sdfs[idat_id][[1]]

    png(filename = paste0(dout, "/qc_plotRedGrnQQ_", idat_id, ".png"), type = "cairo")
    sesameQC_plotRedGrnQQ(sdf)
    dev.off()
    
    png(filename = paste0(dout, "/qc_sesameQC_plotIntensVsBetas_", idat_id, ".png"), type = "cairo")
    sesameQC_plotIntensVsBetas(sdf)
    dev.off()
}

png(filename = paste0(dout, "/qc_sesameQC_plotHeatSNPs.png"), type = "cairo")
sesameQC_plotHeatSNPs(sdfs)  # plots SNP probes and can be used to detect sample swaps
dev.off()

