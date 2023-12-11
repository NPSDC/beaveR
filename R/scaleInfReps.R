# Modified from https://github.com/mikelove/fishpond/blob/master/R/helper.R
#' Compute Size Factors
#'
#' A helper function to compute size factors using median ratio method.
#' First, counts are corrected per row using the effective lengths
#' (for gene counts, the average transcript lengths), then scaled
#' per column to the geometric mean sequence depth, and finally are
#' adjusted per-column up or down by the median ratio size factor to
#' minimize systematic differences across samples.
#'
#' @param tse a TreeSummarizedExperiment with: \code{infReps} a list of
#' inferential replicate count matrices, \code{counts} the
#' estimated counts matrix, and \code{length} the effective
#' lengths matrix
#' @param type whether to use txp or gene for computing size factor(default txp)
#' @param lengthCorrect whether to use effective length correction
#' (default is TRUE)
#' @param meanDepth (optional) user can
#' specify a different mean sequencing depth. By default
#' the geometric mean sequencing depth is computed
#' @param sfFun (optional) size factors function. An
#' alternative to the median ratio can be provided here to adjust
#' the scaledTPM so as to remove remaining library size differences.
#' @param minCount for internal filtering, the minimum count
#' @param minN for internal filtering, the minimum sample size
#' at \code{minCount}
#' @param quiet display no messages
#' @param force boolean to forcefully compute size factors (default FALSE)
#' @return A TreeSummarizedExperiment object where size factors have been
#' computed and added as a column \code{sf} to metadata
#'
#' @examples
#'
#' # path to example data
#' dir <- system.file("extdata/brain_sim_nodtu_small_example", package = "beaveR")
#' # path to file output by TreeTerminus
#' clustFile <- file.path(dir, "cluster_nwk.txt")
#' # path to Salmon quantified files
#' quantDir <- file.path(dir, "out_sal")
#' samples <- as.vector(outer(c(1:6), c(1, 2), function(x, y) paste(x, y, sep = "_")))
#' quantFiles <- file.path(quantDir, samples, "quant.sf")
#' coldata <- data.frame(files = quantFiles, names = samples, condition = factor(rep(1:2, each = 6)))
#' tse <- buildTSE(treeTermFile = clustFile, coldata = coldata)
#' tse <- computeSizeFactors(tse)
#'
#' @export
#' @importFrom TreeSummarizedExperiment rowTree
#' @importFrom tximeta summarizeToGene
computeSizeFactors <- function(tse, type = "txp", lengthCorrect = TRUE,
                               meanDepth = NULL, sfFun = NULL,
                               minCount = 10, minN = 3,
                               quiet = FALSE, force = FALSE) {
  if (!interactive()) {
    quiet <- TRUE
  }

  if (!(is(tse, "TreeSummarizedExperiment"))) {
    stop("tse does not belong to the TreeSummarizedExperiment")
  }

  if (!is.null(metadata(tse)$sf) & !force) {
    if (is.numeric(metadata(tse)$sf)) stop("size factors already present")
  }

  if (!type %in% c("txp", "gene")) {
    stop("type can be only txp or gene")
  }

  metadata(tse)[["type"]] <- type
  infRepIdx <- grep("infRep", assayNames(tse))
  infRepError(infRepIdx)

  tree <- rowTree(tse)
  l <- length(tree$tip.label)

  stopifnot(all(tree$tip.label == rownames(tse)[1:l]))
  tseT <- tse
  tseT <- tse[1:l, ]

  if (type == "gene") {
    gse <- summarizeToGene(tseT)
  }

  counts <- if (type == "txp") assays(tseT)[["counts"]] else assays(gse)[["counts"]]
  infReps <- if (type == "txp") assays(tseT)[infRepIdx] else assays(gse)[infRepIdx]
  length <- if (type == "txp") assays(tseT)[["length"]] else assays(gse)[["length"]]
  nreps <- length(infReps)
  if (is.null(meanDepth) & is.null(sfFun)) {
    meanDepth <- exp(mean(log(colSums(counts))))
  }
  sizeMat <- matrix(nrow = nreps, ncol = ncol(tse))
  if (is.null(length)) {
    if (lengthCorrect) {
      if (!quiet) message("not correcting for feature length (lengthCorrect=FALSE)")
    }
    lengthCorrect <- FALSE
  }
  for (k in seq_len(nreps)) {
    if (!quiet) svMisc::progress(k, max.value = nreps, init = (k == 1), gui = FALSE)
    if (lengthCorrect) {
      # new length bias correction matrix centered on 1
      length <- length / exp(rowMeans(log(length)))
      # a temporary matrix 'cts' which will store
      # the inferential replicate counts
      cts <- infReps[[k]] / length
    } else {
      # for 3' tagged scRNA-seq for example, don't length correct
      cts <- infReps[[k]]
    }
    # if size factors function not provided...
    if (is.null(sfFun)) {
      # divide out the column sum, then set all to the meanDepth

      cts <- t(t(cts) / colSums(cts)) * meanDepth
      # filtering for calculting median ratio size factors
      use <- rowSums(infReps[[k]] >= minCount) >= minN
      loggeomeans <- rowMeans(log(cts[use, ]))
      sf <- apply(cts[use, ], 2, function(s) {
        exp(median((log(s) - loggeomeans)[is.finite(loggeomeans)]))
      })
    } else {
      sf <- sfFun(cts)
    }

    infReps[[k]] <- t(t(cts) / sf)
    sizeMat[k, ] <- sf
  }
  metadata(tse)$sf <- sizeMat
  tse
}

# Modified from https://github.com/mikelove/fishpond/blob/master/R/helper.R
#' Scale inferential replicate counts
#'
#' A helper function to scale the inferential replicates
#' to the mean sequencing depth. The scaling takes into account
#' a robust estimator of size factor (median ratio method is used).
#' First, counts are corrected per row using the effective lengths
#' (for gene counts, the average transcript lengths), then scaled
#' per column to the geometric mean sequence depth, and finally are
#' adjusted per-column up or down by the median ratio size factor to
#' minimize systematic differences across samples.
#'
#' @param tse a TreeSummarizedExperiment with: \code{infReps} a list of
#' inferential replicate count matrices, \code{counts} the
#' estimated counts matrix, and \code{length} the effective
#' lengths matrix
#' @param lengthCorrect whether to use effective length correction
#' (default is TRUE)
#' @param saveMeanScaled store the mean of scaled inferential
#' replicates as an assay 'meanScaled'
#' @param meanDepth (optional) user can
#' specify a different mean sequencing depth. By default
#' the geometric mean sequencing depth is computed
#' @param quiet display no messages
#' @param szMat matrix of size factors computed for each inf replicate (default NULL)
#' @param force boolean to forcefully scale replicates (default FALSE)
#'
#' @return a TreeSummarizedExperiment with the inferential replicates
#' as scaledTPM with library size already corrected (no need for further
#' normalization). A column \code{log10mean} is also added which is the
#' log10 of the mean of scaled counts across all samples and all inferential
#' replicates.
#'
#' @examples
#'
#' # path to example data
#' dir <- system.file("extdata/brain_sim_nodtu_small_example", package = "beaveR")
#' # path to file output by TreeTerminus
#' clustFile <- file.path(dir, "cluster_nwk.txt")
#' # path to Salmon quantified files
#' quantDir <- file.path(dir, "out_sal")
#' samples <- as.vector(outer(c(1:6), c(1, 2), function(x, y) paste(x, y, sep = "_")))
#' quantFiles <- file.path(quantDir, samples, "quant.sf")
#' coldata <- data.frame(files = quantFiles, names = samples, condition = factor(rep(1:2, each = 6)))
#' tse <- buildTSE(treeTermFile = clustFile, coldata = coldata)
#' tse <- computeSizeFactors(tse)
#' tse <- scInfReps(tse)
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' assayNames assayNames<- assay assay<- assays assays<- mcols mcols<-
#' @importFrom TreeSummarizedExperiment rowTree
#' @importFrom stats median
scInfReps <- function(tse, szMat = NULL, lengthCorrect = TRUE,
                         saveMeanScaled = FALSE, quiet = FALSE,
                         force = FALSE, meanDepth = NULL) {
  if (!(is(tse, "SummarizedExperiment"))) {
    stop("tse does not belong to the SummarizedExperiment or the classes that inherit it")
  }

  if (metadata(tse)$infRepsScaled & !force) {
    stop("inferential replicates are already scaled")
  }

  sizeMat <- metadata(tse)$sf
  if (is.null(szMat)) {
    if (is.null(sizeMat)) {
      stop("size factors neither pre-computed nor provided")
    }
  } else {
    if (!is.numeric(szMat)) {
      stop("size factor not numeric")
    }
    sizeMat <- szMat
  }

  infRepIdx <- grep("infRep", assayNames(tse))
  infRepError(infRepIdx)
  infReps <- assays(tse)[infRepIdx]
  nreps <- length(infReps)
  length <- assays(tse)[["length"]]

  if(is.null(metadata(tse)[["txpsAnn"]])) {
      stop("txpsAnn missing in metadata")
  }
  l <- nrow(metadata(tse)[["txpsAnn"]])

  if (!all(dim(sizeMat) == c(nreps, ncol(tse)))) {
    stop("Incorrect dimensions of size factor")
  }

  if (is.null(meanDepth) & is.null(szMat)) {
    meanDepth <- exp(mean(log(colSums(assays(tse)[["counts"]][1:l, ]))))
  }

  metadata(tse)$sf <- sizeMat
  means <- matrix(nrow = nrow(tse), ncol = nreps)

  if (nrow(sizeMat) != nreps | ncol(sizeMat) != ncol(tse)) {
    stop("Invalid dimensions of sizeMat")
  }

  if (is.null(length)) {
    if (lengthCorrect) {
      if (!quiet) message("not correcting for feature length (lengthCorrect=FALSE)")
    }
    lengthCorrect <- FALSE
  }
  if (is.null(szMat)) message("Setting inf rep depth to mean depth")
  for (k in seq_len(nreps)) {
    if (!quiet) svMisc::progress(k, max.value = nreps, init = (k == 1), gui = FALSE)
    if (lengthCorrect) {
      length <- length / exp(rowMeans(log(length)))
      cts <- infReps[[k]] / length
    } else {
      cts <- infReps[[k]]
    }
    if (is.null(szMat)) {
      cts <- t(t(cts) / colSums(cts[1:l, ]) * meanDepth)
    }
    sf <- sizeMat[k, ]
    infReps[[k]] <- t(t(cts) / sf)
    means[, k] <- rowMeans(infReps[[k]])
  }

  if (!quiet) message("")

  assays(tse)[infRepIdx] <- infReps
  mcols(tse)$log10mean <- log10(rowMeans(means) + 1)
  metadata(tse)$infRepsScaled <- TRUE

  if (saveMeanScaled) {
    infRepsArray <- abind::abind(as.list(infReps), along = 3)
    meanScaled <- apply(infRepsArray, 1:2, mean)
    assays(tse)[["meanScaled"]] <- meanScaled
  }
  tse
}

infRepError <- function(infRepIdx) {
  if (length(infRepIdx) == 0) {
    stop("there are no inferential replicates in the assays of 'y'")
  }
}

# #function generator
# defunct = function(msg = "This function is depreciated") function(...) return(stop(msg))

#' @description
#' `r lifecycle::badge("deprecated")`

#' @rdname scInfReps
#' @export
scaleInfReps = function(tse, szMat = NULL, lengthCorrect = TRUE,
                         saveMeanScaled = FALSE, quiet = FALSE,
                         force = FALSE, meanDepth = NULL) {
  lifecycle::deprecate_stop("0.99.1", "scaleInfReps(...)", "scInfReps()")
}
