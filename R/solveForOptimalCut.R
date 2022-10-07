#' Finds a cut and its optimal value on solving the objective function that
#' minimizes(maximizes) the sum for that metric over all nodes in the cut
#'
#' @param tse TreeSumarizedExperiment obtained as the output of running
#' \code{buildTSE}
#' @param metVec numeric vector representing a metric for all nodes in the tree
#' @param type either 'min' or 'max', to minimize or maximize the sum of
#' the metric in the objective function
#'
#' @return List containing the nodes in the cut and the corresponding optimal
#' value
#' @export
#' @examples
#' example(buildTSE)
#' gamma <- 0.1
#' tree <- TreeSummarizedExperiment::rowTree(tse)
#' descSize <- sapply(phangorn::Descendants(tree, seq(nrow(tse))), length)
#' metric <- (mcols(tse)[['meanInfRV']] + ape::node.depth(tree, 2)*gamma) *descSize
#' objS <- solveForOptimalCut(tse, metVec = metric, type = 'min')
#' print(objS[['optVal']])
#' print(length(objS[['cut']]))

#'
#' @importFrom methods is
solveForOptimalCut <- function(tse, metVec, type = "min") {
    if (!is(tse, "TreeSummarizedExperiment")) {
        stop("TSE should be of class experiment")
    }

    if (nrow(tse) != length(metVec)) {
        stop("metVec should be of length equal to number of rows in tse")
    }

    if (!type %in% c("min", "max")) {
        stop("type has to be either min or max")
    }
    findOptSum <- function(tree, metVec, node, type = "min") {
        # print(node)
        if (pkg.globals[["globVec"]][node] != -100) {
            return(pkg.globals[["globVec"]][node])
        }
        if (node <= length(tree$tip)) {
            pkg.globals[["globVec"]][node] <<- metVec[node]
            return(metVec[node])
        }
        children <- phangorn::Descendants(tree, node, "child")
        if (type == "min")
            val <- min(metVec[node], sum(sapply(children, function(child) findOptSum(tree,
                metVec, child, type = type)))) else val <- max(metVec[node], sum(sapply(children, function(child) findOptSum(tree,
            metVec, child, type = type))))
        pkg.globals[["globVec"]][node] <<- val
        return(val)
    }

    findCut <- function(tree, optVals, metVec, node) {
        if (optVals[node] == metVec[node]) {
            return(node)
        }
        children <- phangorn::Descendants(tree, node, "child")
        return(unlist(sapply(children, function(child) findCut(tree, optVals, metVec, child))))
    }

    pkg.globals <- new.env(parent = emptyenv())
    pkg.globals[["globVec"]] <- rep(-100, length(metVec))
    tree <- TreeSummarizedExperiment::rowTree(tse)

    node <- length(tree$tip) + 1
    optVal <- findOptSum(tree, metVec, node, type)
    optVals <- pkg.globals[["globVec"]]

    cut <- findCut(tree, optVals, metVec, length(tree$tip) + 1)
    return(list(cut = cut, optVal = optVal))
}
