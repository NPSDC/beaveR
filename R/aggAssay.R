performMatCheck <- function(seM) {
    if(!(is(seM, "matrix") | is(seM, "dgCMatrix"))) {
        stop("input matrix should be of class dgCMatrix or matrix")
    }
    if(is(seM, "matrix")) {
        if(!is.numeric(seM)) {
            stop("input matrix should be numeric")
        }
    }
}

performIndCheck <- function(seM, inds, type = "row") {
    if(!is(inds, "list")) {
        stop("colInds has to be list")
    }
    if(type=="row") {
        if(any(unlist(inds) > nrow(seM))) {
            stop("Atleast 1 index in rowInds is greater than the total number of rows")
        }
    }
    else {
        if(any(unlist(inds) > ncol(seM))) {
            stop("Atleast 1 index in colInds is greater than the total number of columns")
        }
    }

}

getDimNames <- function(seM, inds, type = "row") {
    d2name <- names(inds)
    if(is.null(names(inds))) {
        d2name <- paste("Group", seq_along(inds), sep="")
    }
    dnames <- list(d2name, colnames(seM))
    if(type!="row")
        dnames <- list(rownames(seM), d2name)
    dnames
}

createMatDf <- function(seM, inds, type="row") {
    nrows <- ifelse(type=="row", length(inds), nrow(seM))
    ncols <- ifelse(type=="col", ncol(seM), length(inds))

    df <- matrix(0, nrow = nrows, ncol = ncols)
    if(is(seM, "dgCMatrix")) {
        df <- Matrix::Matrix(0, nrow = nrows, ncol = ncols, sparse = T)
    }
    df
}

performColMeanAgg <- function(seM, colInds = NULL) {
    performMatCheck(seM)
    if(is.null(colInds)) {
        return(seM)
    }
    performIndCheck(seM, colInds, type= "col")
    df <- createMatDf(seM, colInds, type="col")

    for(i in seq_along(colInds)) {
        if(length(colInds[[i]])!=1) {
            if(is(seM, "dgCMatrix")){
                df[,i] <- Matrix::rowMeans(seM[,colInds[[i]]])
            }
            else {
                df[,i] <- rowMeans(seM[,colInds[[i]]])
            }
        }
        else {
            df[,i] <- seM[,colInds[[i]]]
        }
    }
    dimnames(df) <- getDimNames()
    return(df)
}

performRowSumAgg <- function(seM, rowInds = NULL)
{
    if(!(is(seM, "matrix") | is(seM, "dgCMatrix"))) {
        stop("input matrix should be of class dgCMatrix or matrix")
    }
    if(is(seM, "matrix")) {
        if(!is.numeric(seM)) {
            stop("input matrix should be numeric")
        }
    }
    if(is.null(rowInds)) {
        return(seM)
    }
    if(!is(rowInds, "list")) {
        stop("rowInds has to be list")
    }
    if(any(unlist(rowInds) > nrow(seM))) {
        stop("Atleast 1 index in rowInds is greater than the total rows")
    }
    df <- matrix(0, nrow = length(rowInds), ncol = ncol(seM))
    if(is(seM, "dgCMatrix"))
        df <- Matrix::Matrix(0, nrow = length(rowInds), ncol = ncol(seM), sparse = T)
    for(i in seq_along(rowInds)){
        if(length(rowInds[[i]])!=1) {
            if(is(seM, "dgCMatrix")){
                df[i,] <- Matrix::colSums(seM[rowInds[[i]],])
            }
            else {
                df[i,] <- colSums(seM[rowInds[[i]],])
            }
        }
        else {
            df[i,] <- seM[rowInds[[i]],]
        }
    }
    if(!is.null(names(rowInds)))
        rNames <- names(rowInds)
    else
        rNames <- paste("Group", seq_along(rowInds), sep="")

    rownames(df) <- rNames
    colnames(df) <- colnames(seM)
    return(df)
}

# aggAssay <- function(tree, nodeID, seMat, groupInds = NULL) {
#     if(!is(tree, "phylo"))
#         stop("tree does not belong to class phylo")
#     if(!is.numeric(nodeID))
#     {
#         nodeID <- as.numeric(nodeID)
#         if(sum(is.na(nodeID)) > 0)
#             stop("Node ids contain a non numeric")
#     }
#     if(!(is(seMat, "matrix") | is(seMat, "dgCMatrix"))) {
#         stop("Matrix or sparse matrix needed for seMat")
#     }
#     mat <- matrix(0, nrow=0, ncol=ncol(se_counts))
#     if(is(seMat, "dgCMatrix")) {
#         mat <- Matrix(0, nrow=0, ncol=ncol(se_counts))
#     }
#
#     leaves <- which(nodeID <= nrow(seMat))
#     innNodesInds <- which(nodeID > nrow(seMat))
#
#     lInds <- Descendants(tree, nodeID[innNodesInds], type = "tips")
#     names(lInds) <- as.character(nodeID[innNodesInds])
#     ls <- sapply(lInds, length)
#
#     if(length(leaves) > 0) {
#         mat <- rbind(mat, performRowAgg(seMat[nodeID[leaves],]))
#     }
#
#     if(length(innNodes) > 0) {
#         mat <- rbind(mat, performRowAgg(seMat, lInds))
#     }
#
#     mat <- performColAgg(mat, group_inds)
#
#     return(mat)
# }
#
# computeAggNodes <- function(tree, nodeID, se_counts, group_inds = list(c(1),c(2))) {
#     performColAgg <- function(counts, group_inds) {
#         if(!is.null(dim(counts)))
#             counts <- colSums(counts)
#         vals <- sapply(group_inds, function(x) mean(counts[x]))
#         vals
#     }
#
#     if(!is.numeric(nodeID))
#     {
#         nodeID <- as.numeric(nodeID)
#         if(sum(is.na(nodeID)) > 0)
#             stop("Node ids contain a non numeric")
#     }
#     df <- matrix(0, nrow = length(nodeID), ncol = length(group_inds))
#     leaves <- which(nodeID <= nrow(se_counts))
#     innNodes <- which(nodeID > nrow(se_counts))
#
#     for(i in seq_along(leaves))
#         df[i,] <- performColAgg(se_counts[i,], group_inds)
#
#     lInds <- Descendants(tree, nodeID[innNodes], type = "tips")
#     if(is.null(i))
#         i <- 0
#     for(j in seq_along(lInds))
#         df[i+j,] <- performColAgg(se_counts[lInds[[j]],], group_inds)
#     return(df)
# }
