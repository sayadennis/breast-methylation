library(stringr)
library(sesame)
library(jsonlite)
library(ggplot2)

dout <- "/projects/p30791/methylation/differential_methylation"

hoxa_probes <- list()

for (i in c(1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 13)) {
  hoxa_probes[[paste0("HOXA", i)]] <- names(
    sesameData_getProbesByGene(
      paste0("HOXA", i),
      platform = "EPIC",
      genome = "hg38",
      upstream = 10000,
      downstream = 10000
      # promoter=TRUE
    )
  )
  print(paste0(
    "Finished for HOXA", i, " - ",
    length(hoxa_probes[[paste0("HOXA", i)]]), " probes found!"
  ))
}

json_data <- toJSON(hoxa_probes, pretty = TRUE)

write(json_data, file = paste0(dout, "/hoxa_probes.json"))
