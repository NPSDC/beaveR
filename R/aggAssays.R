aggAssays <- function(tree, se)
{
    if(!is(tree, "phylo")) {
        stop("tree should be of class phylo")
    }
    if(!(is(se, "SummarizedExperiment") | is(se, "SingleCellExperiment"))) {
        stop("se should be either SummarizedExperiment or SingleCellExperiment")
    }
    assays <- SummarizedExperiment::assays
    assayNames <- SummarizedExperiment::assayNames
    rowData <- SummarizedExperiment::rowData
    colData <- SummarizedExperiment::colData
    if(is(se, "SingleCellExperiment")) {
        assays <- SingleCellExperiment::assays
        assayNames <- SingleCellExperiment::assayNames
        rowData <-SingleCellExperiment::rowData
        colData <- SingleCellExperiment::colData
    }

    innNodes <- nrow(se)+1:tree$Nnode
    assaysList <- vector(mode = "list", length(assays(se))-1)
    lInd <- which(assayNames(se) != "length")
    print(length(lInd))
    assaysList <- lapply(assayNames(se)[lInd], function(n) {
        print(n)
        agg <- aggAssay(tree, c(1:nrow(se),innNodes), assays(se)[[n]])
    })
    names(assaysList) <- assayNames(se)[lInd]

    if(is(se, "SummarizedExperiment")) {
        y <- SummarizedExperiment::SummarizedExperiment(assays = assaysList,
                                  colData = colData(se)
                                  )
    } else {
        y <- SingleCellExperiment::SingleCellExperiment(assays = assaysList,
                                                        colData = colData(se)
        )
    }
    metadata <- metadata(se)
    metadata[["txpsAnn"]] = rowData(se)
    metadata(y) <- metadata
    y
}
if(F) {
   a=1
} else {
    a=2
}
