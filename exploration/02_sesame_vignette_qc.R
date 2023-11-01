library(sesame)

#### Calculate Quality Metrics ####

## calculate metrics on all IDATs in a specific folder
qcs = openSesame(idat_dir, prep="", func=sesameQC_calcStats)

sdfs <- sesameDataGet("EPIC.5.SigDF.normal")[1:2] # get two examples
## only compute signal detection stats
qcs = openSesame(sdfs, prep="", func=sesameQC_calcStats, funs="detection")
qcs[[1]]

sesameQC_getStats(qcs[[1]], "frac_dt")

## combine a list of sesameQC into a data frame
head(do.call(rbind, lapply(qcs, as.data.frame)))

sdf <- sesameDataGet('EPIC.1.SigDF')
qc = openSesame(sdf, prep="", func=sesameQC_calcStats, funs=c("detection"))
## equivalent direct call
qc = sesameQC_calcStats(sdf, c("detection"))
qc

#### Rank Quality Metrics ####

sdf <- sesameDataGet('EPIC.1.SigDF')
qc <- sesameQC_calcStats(sdf, "intensity")
qc

sesameQC_rankStats(qc, platform="EPIC")

#### Quality Control Plots ####

## sesameQC_plotBar() takes a list of sesameQC objects and creates bar plot for each metric calculated.

## sesameQC_plotRedGrnQQ() graphs the dye bias between the two color channels.

## sesameQC_plotIntensVsBetas() plots the relationship between β values and signal intensity and can be used to diagnose artificial readout and influence of signal background.

## sesameQC_plotHeatSNPs() plots SNP probes and can be used to detect sample swaps.

