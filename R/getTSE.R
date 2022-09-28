#' @param txps
#'
#' @param y
#' @param ...
#'
#' @importFrom fishpond labelKeep
getTxps <- function(txps, y, ...) {
    nMissParams <- c(!missing(txps),!missing(y))
    if (!all(nMissParams)) {
        stop("Missing Parameters to the function")
    }
    if (is.null(txps)) {
        y <- labelKeep(y, ...)
        if (sum(mcols(y)[["keep"]]) == 0) {
            stop("No txps match the filtering criteria, change either minN or minCount")
        }
        txps <- rownames(y)[mcols(y)[["keep"]]]
    }
    else {
        if (!all(txps %in% rownames(y)))
            stop("txps provided missing in the data")
    }
    return(txps)
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
    txps <- getTxps(txps, se, ...)
    treeMerged <- mergeTrees(treeTerm, tnames = rownames(se))
    mb <- mergeLeaves(treeCons, se[txps,])

}
