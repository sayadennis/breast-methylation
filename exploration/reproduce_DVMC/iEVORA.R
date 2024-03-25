#### iEVORA.R
#### Author: Andrew E Teschendorff (a.teschendorff@ucl.ac.uk)
#### Date: 26th Nov.2015
#### Copyright 2015 Andrew Teschendorff

library(qvalue)
library(stats)

doDV <- function(tmp.v, pheno.v) {
  co.idx <- which(pheno.v == 0)
  ca.idx <- which(pheno.v == 1)
  bt.o <- bartlett.test(x = tmp.v, g = pheno.v, na.action = na.omit)
  pv <- bt.o$p.value
  logR <- log2(var(tmp.v[ca.idx], na.rm = TRUE) / var(tmp.v[co.idx], na.rm = TRUE))
  avCA <- mean(tmp.v[ca.idx], na.rm = TRUE)
  avCO <- mean(tmp.v[co.idx], na.rm = TRUE)
  out.v <- c(logR, pv, avCA, avCO)
  names(out.v) <- c("log(V1/V0)", "P(BT)", "Av1", "Av0")
  return(out.v)
}


doTT <- function(tmp.v, pheno.v) {
  tt.o <- t.test(tmp.v ~ pheno.v)
  out.v <- c(-tt.o$stat, tt.o$p.val)
  names(out.v) <- c("t", "P")
  return(out.v)
}

iEVORA <- function(data.m, pheno.v, thDV = 0.001, thDM = 0.05, return_top = FALSE) {
  statDVC.m <- t(apply(data.m, 1, doDV, pheno.v))
  print("Estimated DV statistics")
  qvDVC.v <- qvalue(statDVC.m[, 2])$qval
  statDMC.m <- t(apply(data.m, 1, doTT, pheno.v))
  print("Preparing output")
  if (return_top == FALSE) {
    DVMC.m <- cbind(
      statDMC.m,
      statDVC.m[, c(3:4, 1:2)],
      qvDVC.v
    )
    colnames(DVMC.m) <- c("t", "P(TT)", "Av1", "Av0", "log[V1/V0]", "P(BT)", "q(BT)")
    rownames(DVMC.m) <- rownames(statDMC.m)
  } else if ((ntop > 0) && (return_top == TRUE)) {
    DVMC.m <- cbind(
      statDMC.m[tmp.s$ix[1:ntop], ],
      statDVC.m[dvc.idx[tmp.s$ix[1:ntop]], c(3:4, 1:2)],
      qvDVC.v[dvc.idx[tmp.s$ix[1:ntop]]]
    )
    colnames(DVMC.m) <- c("t", "P(TT)", "Av1", "Av0", "log[V1/V0]", "P(BT)", "q(BT)")
    rownames(DVMC.m) <- rownames(statDMC.m)[tmp.s$ix[1:ntop]]
  } else {
    print("No DVMCs in this data.")
  }

  return(DVMC.m)
}
