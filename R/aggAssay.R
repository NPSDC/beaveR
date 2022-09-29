aggAssay <- function(tree, nodeIDs, seMat, groupInds = NULL) {
    if (!is(tree, "phylo"))
        stop("tree does not belong to class phylo")
    if (!is.numeric(nodeIDs))
    {
        nodeIDs <- suppressWarnings(as.numeric(nodeIDs))
        if (any(is.na(nodeIDs)))
            stop("Node ids contain a non numeric")
    }
    if (any(nodeIDs > length(tree$tip) + tree$Nnode))
        stop("some nodeIDs are invalid")
    performMatCheck(seMat)
    mat <- createMat(seMat, list(), type = "row")

    leaves <- which(nodeIDs <= nrow(seMat))
    innNodesInds <- which(nodeIDs > nrow(seMat))

    if (length(leaves) > 0) {
        mat <- rbind(mat, seMat[nodeIDs[leaves], ])
        rownames(mat) <- rownames(seMat)[nodeIDs[leaves]]
    }

    if (length(innNodesInds) > 0) {
        lInds <-
            phangorn::Descendants(tree, nodeIDs[innNodesInds], type = "tips")
        names(lInds) <-
            paste("Node", as.character(nodeIDs[innNodesInds]), sep = "")
        mat <- rbind(mat, performRowSumAgg(seMat, lInds))
    }
    mat <- performColMeanAgg(mat, groupInds)
    mat
}

performColMeanAgg <- function(seM, colInds = NULL) {
    performMatCheck(seM)
    if (is.null(colInds)) {
        return(seM)
    }
    performIndCheck(seM, colInds, type = "col")
    mat <- createMat(seM, colInds, type = "col")

    for (i in seq_along(colInds)) {
        if (length(colInds[[i]]) != 1) {
            if (is(seM, "dgCMatrix")) {
                mat[, i] <- Matrix::rowMeans(seM[, colInds[[i]]])
            }
            else {
                mat[, i] <- rowMeans(seM[, colInds[[i]]])
            }
        }
        else {
            mat[, i] <- seM[, colInds[[i]]]
        }
    }
    dimnames(mat) <- getDimNames(seM, colInds, type = "col")
    mat
}

performRowSumAgg <- function(seM, rowInds = NULL)
{
    performMatCheck(seM)
    if (is.null(rowInds)) {
        return(seM)
    }
    performIndCheck(seM, rowInds, type = "row")
    mat <- createMat(seM, rowInds, type = "row")

    for (i in seq_along(rowInds)) {
        if (length(rowInds[[i]]) != 1) {
            if (is(seM, "dgCMatrix")) {
                mat[i, ] <- Matrix::colSums(seM[rowInds[[i]], ])
            }
            else {
                mat[i, ] <- colSums(seM[rowInds[[i]], ])
            }
        }
        else {
            mat[i, ] <- seM[rowInds[[i]], ]
        }
    }
    dimnames(mat) <- getDimNames(seM, rowInds, type = "row")
    mat
}

performMatCheck <- function(seM) {
    if (!(is(seM, "matrix") | is(seM, "dgCMatrix"))) {
        stop("input matrix should be of class dgCMatrix or matrix")
    }
    if (is(seM, "matrix")) {
        if (!is.numeric(seM)) {
            stop("input matrix should be numeric")
        }
    }
}

performIndCheck <- function(seM, inds, type = "row") {
    if (!is(inds, "list")) {
        stop("colInds has to be list")
    }
    if (type == "row") {
        if (any(unlist(inds) > nrow(seM))) {
            stop("Atleast 1 index in rowInds is greater than the total number of rows")
        }
    }
    else {
        if (any(unlist(inds) > ncol(seM))) {
            stop("Atleast 1 index in colInds is greater than the total number of columns")
        }
    }

}

getDimNames <- function(seM, inds, type = "row") {
    d2name <- names(inds)
    if (is.null(names(inds))) {
        d2name <- paste("Group", seq_along(inds), sep = "")
    }
    dnames <- list(d2name, colnames(seM))
    if (type != "row")
        dnames <- list(rownames(seM), d2name)
    dnames
}

createMat <- function(seM, inds, type = "row") {
    nrows <- ifelse(type == "row", length(inds), nrow(seM))
    ncols <- ifelse(type == "col", length(inds), ncol(seM))

    mat <- matrix(0, nrow = nrows, ncol = ncols)
    if (is(seM, "dgCMatrix")) {
        mat <- Matrix::Matrix(0,
                              nrow = nrows,
                              ncol = ncols,
                              sparse = T)
    }
    mat
}
