#' @importFrom SummarizedExperiment assays assayNames colData rowData metadata
aggAssays <- function(tree, se)
{
    if(!is(tree, "phylo")) {
        stop("tree should be of class phylo")
    }
    if(!(is(se, "SummarizedExperiment") | is(se, "SingleCellExperiment"))) {
        stop("se should be either SummarizedExperiment or SingleCellExperiment")
    }
    innNodes <- nrow(se)+1:tree$Nnode
    assaysList <- vector(mode = "list", length(assays(se))-1)
    lInd <- which(assayNames(se) != "length")
    print(length(lInd))
    assaysList <- lapply(assayNames(se)[lInd], function(n) {
        print(n)
        agg <- aggAssay(tree, c(1:nrow(se),innNodes), assays(se)[[n]])
        print(dim(agg))
        agg
    })
    names(assaysList) <- assayNames(se)[lInd]

    if(is(se, "SummarizedExperiment")) {
        y <- SummarizedExperiment::SummarizedExperiment(assays = assaysList,
                                  colData = colData(se),
                                  metadata = metadata(se)
                                  )
        metadata(y)[["txpsAnn"]] = rowData(se)
        return(y)
    }
    y <- SingleCellExperiment::SingleCellExperiment(assays = assaysList,
                                                    colData = colData(se),
                                                    metadata = metadata(se)
                                                    )
    metadata(y)[["txpsAnn"]] = rowData(se)
    y
}
