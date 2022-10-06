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

makeTSEFromSE <- function(tree, se) {
    if (!(is(se, "SummarizedExperiment") |
          is(se, "SingleCellExperiment"))) {
        stop("se can only come from SummarizedExperiment/SingleCellExperiment")
    }
    if (!is(tree, "phylo")) {
        stop("tree should be from phylo class")
    }
    l <- length(tree$tip.label)
    if(l != nrow(se)) {
        stop("number of leaves should be equal to number of rows in se")
    }
    if(!all(tree$tip.label %in% rownames(se))) {
        stop("leaf set of tree is not same as the transcripts in se")
    }
    se <- se[tree$tip.label,]
    seAgg <- aggAssays(tree, se)
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
        assays = assays(seAgg),
        colData = colData(seAgg),
        metadata = metadata(seAgg),
        rowTree = tree
    )
    mcols(tse) <- mcols(seAgg)
    tse
}

#' @param treeTermFile
#'
#' @param coldata
#' @param txps
#' @param ...
#'
#' @export
makeTSE <- function(treeTermFile,
                   coldata,
                   txps = NULL,
                   ...) {
    if (!file.exists(treeTermFile)) {
        stop(paste("the file", treeTermFile, "does not exist"))
    }
    message("reading tree")
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

    tse <- makeTSEFromSE(treeMerged, se)
    tse <- fishpond::computeInfRV(tse, meanVariance = FALSE)
    tse
}
