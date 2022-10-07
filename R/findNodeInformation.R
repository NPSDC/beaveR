#' @importFrom phangorn Descendants
findGenes <- function(tree, node, txpsAnn) {
    desc <- unlist(Descendants(tree, node, "tips"))
    genes <- txpsAnn[tree$tip.label[desc], "gene_id"]
    genes
}

#' @export
#' @importFrom phangorn Descendants
#' @importFrom methods is
#' @importFrom S4Vectors metadata
findNodeInformation <- function(tse, node, type = NULL) {
    if (!is(tse, "TreeSummarizedExperiment")) {
        stop("tse is not Tree Summarized Object")
    }
    if (!is.numeric(node)) {
        stop("node has to be numeric")
    }
    if (node > nrow(tse)) {
        stop("node index is invalid since it is larger than the number of rows in tse ")
    }

    tree <- TreeSummarizedExperiment::rowTree(tse)
    df <-
        data.frame(nodeInd = node, meanInfRV = SummarizedExperiment::mcols(tse)[["meanInfRV"]][node])
    txpsAnn <- metadata(tse)[["txpsAnn"]]


    if ("gene_id" %in% colnames(txpsAnn)) {
        genes <- list(findGenes(tree, node, txpsAnn))
        df <- cbind(df, genes = I(genes))
    }

    if(!is.null(type)) {
        if (!(type %in% c("tips", "children", "all"))) {
            stop("type can only be - tips, children, all")
        }
        desc <- unlist(Descendants(tree, node, type = type))
        if (length(desc) > 0) {
            desc <- setdiff(desc, node)
            mirv <- SummarizedExperiment::mcols(tse)[["meanInfRV"]][desc]
            dfD <- data.frame(nodeInd = desc, meanInfRV = mirv)

            txpsAnn <- metadata(tse)[["txpsAnn"]]
            if ("gene_id" %in% colnames(txpsAnn)) {
                genes <- lapply(desc, function(n) {
                    findGenes(tree, n, txpsAnn)
                })
                dfD <- cbind(dfD, genes = I(genes))
            }
            df <- rbind(df, dfD)
        }
    }
    nodeType <- rep("Leaf", nrow(df))
    nodeType[df$nodeInd > length(tree$tip.label)]="InnerNode"
    df <- cbind(df, nodeType=nodeType)
    df
}
#
