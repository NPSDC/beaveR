mergeTrees <- function(trees, updateInd = T, tnames = NULL) {
    if (updateInd) {
        ### Because of Salmon is 0 index and R is 1 index
        trees <- lapply(trees, function(tree) {
            tree$tip.label = as.character(as.numeric(tree$tip.label) + 1)
            tree
        })
    }

    nwk <- paste("(", paste(sapply(trees, function(tree) gsub(";", "", ape::write.tree(tree))),
        collapse = ","), ");", sep = "")
    tree <- ape::read.tree(text = nwk)
    if (!is.null(tnames)) {
        if (length(tnames) < length(tree$tip.label)) {
            stop("length of tnames is smaller")
        }
        if (any(as.numeric(tree$tip.label) - 1 > length(tnames))) {
            stop("some tree tips are greater than total number of rows of data frame")
        }
        tree$tip.label <- tnames[as.numeric(tree$tip.label)]
    }
    return(tree)
}
