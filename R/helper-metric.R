# https://github.com/mikelove/fishpond/blob/c1e968f810f2caba2b333a3891d2c3c05e73001b/R/swish.R#L338
#' @importFrom matrixStats rowMedians
getLog2FC <- function(infRepsArray, condition, pc = 5) {
  dims <- dim(infRepsArray)
  cond1 <- condition == levels(condition)[1]
  cond2 <- condition == levels(condition)[2]
  diffs <- matrix(nrow = dims[1], ncol = dims[3])
  for (k in seq_len(dims[3])) {
    diffs[, k] <-
      log2(rowMeans(infRepsArray[, cond2, k]) + pc) - log2(rowMeans(infRepsArray[
        ,
        cond1, k
      ]) + pc)
  }

  # median over inferential replicates
  rowMedians(diffs)
}


#' Computes log fold change on the scaled inferential replicate counts.
#' Scaled count for a given transcript in an inferential replicate is given by
#' dividing the original count by the sum of counts of all transcripts in that
#' inferential replicate
#'
#' @param tse TreeSumarizedExperiment obtained as the output of running
#' \code{buildTSE}
#' @param x the name of the condition variable. A factor with two levels for a
#' two group analysis. The log fold change is computed as non-reference level
#' over reference level
#' @param pc integer pseudocount for computing log fold change
#'
#' @return numeric vector of length equal to number of rows in tse
#'
#' @export
#' @examples
#' example(buildTSE)
#' lfc <- getScaledLFC(tse, "condition")
#' @importFrom SummarizedExperiment colData assays assayNames
#' @importFrom TreeSummarizedExperiment rowTree
#' @importFrom methods is
getScaledLFC <- function(tse, x, pc = 5) {
  if (!(x %in% colnames(colData(tse)))) {
    stop("x is not a valid column")
  }
  if (!(is(tse, "TreeSummarizedExperiment"))) {
    stop("input se object not from the right class")
  }
  cond <- colData(tse)[[x]]
  if (!is(cond, "factor")) {
    stop("x is not factor")
  }

  infRepInds <- grep("infRep", assayNames(tse))
  if (sum(infRepInds) == 0) {
    stop("infReps not present in the assays")
  }
  tree <- rowTree(tse)
  l <- length(tree$tip)
  infReps <- assays(tse)[infRepInds]
  infReps <- abind::abind(as.list(infReps), along = 3)
  mSf <- 0
  for (j in seq(dim(infReps)[3])) {
    sf <- colSums(infReps[1:l, , j])
    mSf <- mean(sf) + mSf
    infReps[, , j] <- t(t(infReps[, , j]) / sf)
  }
  mSf <- mSf / dim(infReps)[3]
  lfc <- getLog2FC(infReps, cond, pc = 5 / mSf)
  lfc
}
