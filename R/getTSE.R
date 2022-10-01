#' @param txps
#'
#' @param y
#' @param ...
#'
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
        txps <-
            rownames(y)[SummarizedExperiment::mcols(y)[["keep"]]]
    }
    else {
        if (!all(txps %in% rownames(y)))
            stop("txps provided missing in the data")
    }
    return(txps)
}

#' @export
makeTSEFromSE <- function(se, tree) {
    if (!(is(se, "SummarizedExperiment") |
          is(se, "SingleCellExperiment"))) {
        stop("se can only come from SummarizedExperiment/SingleCellExperiment")
    }
    if (!is(tree, "phylo")) {
        stop("tree should be from phylo class")
    }
    assays <- SummarizedExperiment::assays
    colData <- SummarizedExperiment::colData
    mcols <- SummarizedExperiment::mcols
    if (is(se, "SingleCellExperiment")) {
        assays <- SingleCellExperiment::assays
        colData <- SingleCellExperiment::colData
    }
    tree$node.label <-
        as.character(paste("Node", length(tree$tip) + 1:tree$Nnode, sep = ""))
    tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
        assays = assays(se),
        colData = colData(se),
        metadata = metadata(se),
        rowTree = tree
    )
    mcols(tse) <- mcols(se)
    tse
}

#' @param treeTermFile
#'
#' @param coldata
#' @param txps
#' @param ...
#'
#' @export
getTSE <- function(treeTermFile,
                   coldata,
                   txps = NULL,
                   ...) {
    if (!file.exists(treeTermFile)) {
        stop(paste("the file", treeTermFile, "does not exist"))
    }
    print("reading tree")
    treeTerm <- ape::read.tree(treeTermFile) ## Reading tree
    if (is.null(treeTerm)) {
        stop(
            paste(
                treeTermFile,
                "could not be parsed, make sure the file follows the correct format expected by ape::read.tree"
            )
        )
    }
    se <- tximeta::tximeta(coldata)
    txpsFilt <- getTxps(txps, se, ...)
    treeMerged <- mergeTrees(treeTerm, tnames = rownames(se))
    treeMerged <- mergeTreeWithSE(treeMerged, txpsFilt)
    se <- se[treeMerged$tip.label, ]

    seAgg <- aggAssays(treeMerged, se)
    # seAgg <- fishpond::computeInfRV(seAgg, meanVariance = FALSE)
    tse <- makeTSEFromSE(seAgg, treeMerged)
    tse <- fishpond::computeInfRV(tse, meanVariance = FALSE)
    tse
}
