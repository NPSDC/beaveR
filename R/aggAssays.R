## https://github.com/mikelove/tximport/blob/092f0a3bf2a1cd2430b7e1c2fa437426e746ea4e/R/helper.R
replaceMissingLength <- function(lengthMat, aveLengthSampGroup) {
  nanRows <- which(apply(lengthMat, 1, function(row) any(is.nan(row))))
  if (length(nanRows) > 0) {
    for (i in nanRows) {
      if (all(is.nan(lengthMat[i, ]))) {
        # if all samples have 0 abundances for all tx, use the simple average
        lengthMat[i, ] <- aveLengthSampGroup[i]
      } else {
        # otherwise use the geometric mean of the lengths from the other samples
        idx <- is.nan(lengthMat[i, ])
        lengthMat[i, idx] <- exp(mean(log(lengthMat[i, !idx]), na.rm = TRUE))
      }
    }
  }
  lengthMat
}


#' @importFrom methods is
#' @importFrom S4Vectors metadata metadata<-
aggAssays <- function(tree, se, groupInds = NULL) {
  if (!is(tree, "phylo")) {
    stop("tree should be of class phylo")
  }
  if (!(is(se, "SummarizedExperiment") |
    is(se, "SingleCellExperiment"))) {
    stop("se should be either SummarizedExperiment or SingleCellExperiment")
  }
  if (!is.null(groupInds) & !is(groupInds, "list")) {
    stop("groupInds is not list")
  }
  assays <- SummarizedExperiment::assays
  assayNames <- SummarizedExperiment::assayNames
  rowData <- SummarizedExperiment::rowData
  colData <- SummarizedExperiment::colData
  l <- length(tree$tip)
  # if (is(se, "SingleCellExperiment")) {
  #     assays <- SingleCellExperiment::assays
  #     assayNames <- SingleCellExperiment::assayNames
  #     rowData <- SingleCellExperiment::rowData
  #     colData <- SingleCellExperiment::colData
  # }

  innNodes <- nrow(se) + 1:tree$Nnode
  reqAssayNames <-
    intersect(
      c("counts", "abundance", assayNames(se)[grep("infRep", assayNames(se))]),
      assayNames(se)
    )

  assaysList <- vector(mode = "list", length(reqAssayNames))
  lInd <- which(assayNames(se) != "length")
  message("Aggregation Started")
  assaysList <- lapply(reqAssayNames, function(n) {
    agg <-
      aggAssay(tree, c(1:nrow(se), innNodes), assays(se)[[n]], groupInds = groupInds)
  })
  names(assaysList) <- reqAssayNames

  # aggregated length
  # Adapted from SummarizeToGene
  # https://github.com/mikelove/tximport/blob/master/R/summarizeToGene.R

  # the next lines calculate a weighted average of transcript length,
  # weighting by transcript abundance.
  # this can be used as an offset / normalization factor which removes length bias
  # for the differential analysis of estimated counts summarized at the group level.
  weightedLength <- aggAssay(tree, innNodes, assays(se)[["abundance"]] * assays(se)[["length"]], groupInds = groupInds)
  lengthMat <- weightedLength / assaysList[["abundance"]][innNodes,]

  # pre-calculate a simple average transcript length
  # for the case the abundances are all zero for all samples.
  # first, average the tx lengths over samples
  aveLengthSamp <- rowMeans(assays(se)[["length"]])

  # then simple average of lengths within groups (not weighted by abundance)
  desc <- phangorn::Descendants(tree, innNodes)
  aveLengthSampGroup <- sapply(desc, function(inds) mean(aveLengthSamp[inds]))

 lengthMat <- replaceMissingLength(lengthMat, aveLengthSampGroup)
 lengthMat <- rbind(assays(se)[["length"]], lengthMat)

  l <- length(assaysList)
  assaysList[[l + 1]] <- lengthMat

  names(assaysList)[l + 1] <- "length"
  assaysList <- assaysList[c(1, 2, l + 1, 3:l)]

  message("Aggregation Ended")
  if (is(se, "SummarizedExperiment")) {
    y <-
      SummarizedExperiment::SummarizedExperiment(assays = assaysList, colData = colData(se))
  } else {
    y <-
      SingleCellExperiment::SingleCellExperiment(assays = assaysList, colData = colData(se))
  }
  metadata <- metadata(se)
  metadata[["txpsAnn"]] <- rowData(se)
  # metadata[["txpsLength"]] <- assays(se)[["length"]]
  metadata(y) <- metadata
  y
}
