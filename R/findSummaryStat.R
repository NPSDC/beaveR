

#' #' @importFrom phangorn Descendants
#' findSummaryStat(tse, node, type = "tips") {
#'     if(!is(tse, "TreeSummarizedExperiment")) {
#'         stop("tse is not Tree Summarized Object")
#'     }
#'     if(!is.numeric(node)) {
#'         stop("node has to be numeric")
#'     }
#'     if(node > nrow(tse)) {
#'         stop("node index is invalid since it is larger than the number of rows in tse ")
#'     }
#'     if(!(type %in% c("tips", "children", "all"))) {
#'         stop("type can only be - tips, children, all")
#'     }
#'     tree <- TreeSummarizedExperiment::rowTree(tse)
#'
#' }
#
