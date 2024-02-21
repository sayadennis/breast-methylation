library(IRanges)
library(sesame)

# nolint start: cyclocomp_linter

sesameData_annoProbes_custom <- function(Probe_IDs, regs = NULL,
                                         collapse = TRUE, chooseOne = FALSE,
                                         column = NULL, sep = ",",
                                         return_ov_probes = FALSE, return_ov_features = FALSE,
                                         out_name = NULL, platform = NULL,
                                         genome = NULL, silent = FALSE) {
  stopifnot(is.character(Probe_IDs))
  if (is.null(platform)) {
    platform <- inferPlatformFromProbeIDs(Probe_IDs, silent = silent)
  }

  if (is.null(regs)) {
    # default to annotate genes
    if (is.null(genome)) {
      genome <- sesameData_check_genome(NULL, platform)
    }
    regs <- sesameData_getTxnGRanges(genome)
    if (is.null(column)) column <- "gene_name"
  }

  gr <- sesameData_getManifestGRanges(platform, genome = genome)
  in_mft <- Probe_IDs %in% names(gr)
  if (sum(!in_mft) > 0) {
    warning(sprintf(
      "%d probes out of manifest were excluded.", sum(!in_mft)
    ))
  }
  probes <- gr[Probe_IDs[in_mft]]
  if (length(probes) == 0) {
    return(probes)
  } # empty GRanges
  hits <- findOverlaps(probes, regs, ignore.strand = TRUE)

  if (return_ov_probes) {
    return(probes[unique(queryHits(hits))])
  } else if (return_ov_features) {
    return(regs[unique(subjectHits(hits))])
  }

  if (is.null(column)) {
    label <- names(regs[subjectHits(hits)])
    if (is.null(out_name)) {
      out_name <- "anno"
    }
  } else {
    stopifnot(column %in% colnames(mcols(regs)))
    label <- mcols(regs[subjectHits(hits)])[[column]]
    if (is.null(out_name)) {
      out_name <- column
    }
  }
  if (collapse) {
    if (chooseOne) {
      pid2label <- vapply(
        split(label, queryHits(hits)),
        function(x) x[1], character(1)
      )
    } else {
      pid2label <- vapply(
        split(label, queryHits(hits)),
        function(x) paste0(unique(x), collapse = sep), character(1)
      )
    }
    mcols(probes)[[out_name]] <- NA
    mcols(probes[as.integer(names(pid2label))])[[out_name]] <- pid2label
  } else {
    unfound <- probes[!(seq_along(probes) %in% queryHits(hits))]
    if (length(unfound) > 0) {
      mcols(unfound)[[out_name]] <- NA
    }
    found <- probes[queryHits(hits)]
    if (length(found) > 0) {
      mcols(found)[[out_name]] <- label
    }
    probes <- sort(c(unfound, found))
  }
  probes
}
# nolint end
