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
    if (is(se, "SingleCellExperiment")) {
        assays <- SingleCellExperiment::assays
        assayNames <- SingleCellExperiment::assayNames
        rowData <- SingleCellExperiment::rowData
        colData <- SingleCellExperiment::colData
    }

    innNodes <- nrow(se) + 1:tree$Nnode
    reqAssayNames <-
        intersect(c("counts", "abundance", assayNames(se)[grep("infRep", assayNames(se))]),
                  assayNames(se))

    assaysList <- vector(mode = "list", length(reqAssayNames))
    lInd <- which(assayNames(se) != "length")
    message("Aggregation Started")
    assaysList <- lapply(reqAssayNames, function(n) {
        agg <-
            aggAssay(tree, c(1:nrow(se), innNodes), assays(se)[[n]], groupInds = groupInds)
    })
    names(assaysList) <- reqAssayNames

    lAssay <- matrix(NA, nrow = dim(assaysList[[1]])[1],
                     ncol = dim(assaysList[[1]])[2])
    lAssay[1:nrow(se),1:ncol(se)] = assays(se)[["length"]]
    l <- length(assaysList)
    assaysList[[l+1]] <- lAssay

    names(assaysList)[l+1] <- "length"
    assaysList <- assaysList[c(1,2,l+1,3:l)]

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
    metadata[["txpsLength"]] <- assays(se)[["length"]]
    metadata(y) <- metadata
    y
}
