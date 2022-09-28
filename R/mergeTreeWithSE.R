mergeTreeWithSE <- function(tree, se) {
    skip("skip")
    txps <- suppressWarnings(as.numeric(tree$tip))
    if(!all(is.na(txps))){
        stop("tree tips contain numeric transcripts")
    }
    missing_txps <- setdiff(rownames(se), tree$tip.label)
    if(length(missing_txps) > 0)
    {
        print(paste("Missing txps", length(missing_txps)))
        remLeaves <- paste(as.character(missing_txps), collapse = ",")
        nwk <- ape::write.tree(tree)
        nwk <- substr(nwk, 2, nchar(nwk)-2)
        nwk <- paste("(", remLeaves, ",", nwk, ");", sep = "")
        tree <- ape::read.tree(text = nwk)
    }
    return(tree)
}
