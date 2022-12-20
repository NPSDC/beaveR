#' @importFrom fishpond labelKeep
getTxps <- function(txps, y, ...) {
  nMissParams <- c(!missing(txps), !missing(y))
  if (!all(nMissParams)) {
    stop("Missing Parameters to the function")
  }
  if (is.null(txps)) {
    y <- labelKeep(y, ...)
    if (sum(SummarizedExperiment::mcols(y)[["keep"]]) == 0) {
      stop("No txps match the filtering criteria, change either minN or minCount")
    }
    txps <- rownames(y)[SummarizedExperiment::mcols(y)[["keep"]]]
  } else {
    if (!all(txps %in% rownames(y))) {
      stop("txps provided missing in the data")
    }
  }
  return(txps)
}

#' @importFrom methods is
buildTSEFromSE <- function(tree, se) {
  if (!(is(se, "SummarizedExperiment") | is(se, "SingleCellExperiment"))) {
    stop("se can only come from SummarizedExperiment/SingleCellExperiment")
  }
  if (!is(tree, "phylo")) {
    stop("tree should be from phylo class")
  }
  l <- length(tree$tip.label)
  if (l != nrow(se)) {
    stop("number of leaves should be equal to number of rows in se")
  }
  if (!all(tree$tip.label %in% rownames(se))) {
    stop("leaf set of tree is not same as the transcripts in se")
  }
  se <- se[tree$tip.label, ]
  seAgg <- aggAssays(tree, se)
  assays <- SummarizedExperiment::assays
  colData <- SummarizedExperiment::colData
  mcols <- SummarizedExperiment::mcols
  # if (is(se, "SingleCellExperiment")) {
  #     assays <- SingleCellExperiment::assays
  #     colData <- SingleCellExperiment::colData
  # }
  tree$node.label <- as.character(paste("Node", length(tree$tip) + 1:tree$Nnode, sep = ""))
  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = assays(seAgg), colData = colData(seAgg),
    metadata = metadata(seAgg), rowTree = tree
  )
  mcols(tse) <- mcols(seAgg)
  metadata(tse)[["infRepsScaled"]] <- FALSE
  tse
}

#' Builds a TreeSummarizedExperiment (TSE) object from salmon quantified files and
#' file of trees given by TreeTerminus. Internally it reads the salmon
#' quantified files using tximeta. The \code{rowTree} is the tree on a transcript set.
#' The transcript set is computed by taking a union of all transcripts that are
#' covered by the trees from TreeTerminus and transcripts that pass the filtering
#' criteria determined using \code{labelKeep} function in \code{fishpond}. The
#' rows consists of all the nodes including the leaf and inner nodes of the
#' \code{rowTree}. The assays consist of (counts, abundance and infRep matrices)
#' with the inner node value of an individual assay representing an aggregation
#' of the descendant leaves for that node in that assay. It requires that Salmon
#' should have been run with either bootstrap or gibbs sampling.
#'
#' @param treeTermFile path to file of trees given by TreeTerminus
#' @param coldata data.frame that is given as an input to tximeta
#' @param txps (Optional) A character vector containing transcripts to which
#' TreeSummarizedExperiment object will be restricted.
#' @param ... arguments passed to \code{labelKeep} function in \code{fishpond}
#'
#' @return TreeSummarizedExperiment Object
#'
#' @examples
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
#'

#' @export
#' @importFrom methods is
#' @importFrom SummarizedExperiment mcols<-
buildTSE <- function(treeTermFile, coldata, txps = NULL, ...) {
  if (!file.exists(treeTermFile)) {
    stop(paste("the file", treeTermFile, "does not exist"))
  }
  message("reading tree")
  treeTerm <- ape::read.tree(treeTermFile) ## Reading tree
  if (is.null(treeTerm)) {
    stop(paste(treeTermFile, "could not be parsed, make sure the file follows the correct format expected by ape::read.tree"))
  }
  se <- tximeta::tximeta(coldata)
  txpsFilt <- getTxps(txps, se, ...)
  treeMerged <- mergeTrees(treeTerm, tnames = rownames(se))
  treeMerged <- mergeTreeWithSE(treeMerged, txpsFilt)
  se <- se[treeMerged$tip.label, ]

  tse <- buildTSEFromSE(treeMerged, se)
  tse <- fishpond::computeInfRV(tse, meanVariance = FALSE)
  tse
}
