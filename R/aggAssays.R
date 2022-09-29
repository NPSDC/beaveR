aggAssays <- function(tree, se)
{
    if(!is(tree, "phylo")) {
        stop("tree should be of class phylo")
    }
    if(!(is(se, "SummarizedExperiment") | is(se, "SingleCellExperiment"))) {
        stop("se should be either SummarizedExperiment or SingleCellExperiment")
    }
    innNodes <- nrow(se)+1:tree$Nnode
    assaysList <- vector(mode = "list", length(assays(se)))
    names(assaysList) <- assayNames(se)
    for(n in names(asList)) {
        if(n == "length") {
            assaysList[[n]] <- assays(se)[[n]]
        }
        else{
            assaysList[[n]] <- aggAssay(tree, c(1:nrow(se),innNodes), assays(se)[[n]])
        }
    }

    if(is(se, "SummarizedExperiment")) {
        y <- SummarizedExperiment::SummarizedExperiment(assays = assaysList,
                                  colData = colData(se),
                                  metadata = metadata(se),
                                  rowData = rowData(se))
        return(y)
    }
    y <- SingleCellExperiment::SingleCellExperiment(assays = assaysList,
                                                    colData = colData(se),
                                                    metadata = metadata(se),
                                                    rowData = rowData(se))
    return(y)

    metadata(y)$infRepsScaled=F
    y
}
