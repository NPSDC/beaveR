mergeTreeWithSE <- function(tree, txpsFilt) {
    txps <- suppressWarnings(as.numeric(tree$tip))
    if (!all(is.na(txps))) {
        stop("tree tips contain numeric transcripts")
    }
    missingTxps <- setdiff(txpsFilt, tree$tip.label)
    if (length(missingTxps) > 0) {
        message(paste("Missing txps", length(missingTxps)))
        remLeaves <- paste(as.character(missingTxps), collapse = ",")
        nwk <- ape::write.tree(tree)
        nwk <- substr(nwk, 2, nchar(nwk) - 2)
        nwk <- paste("(", remLeaves, ",", nwk, ");", sep = "")
        tree <- ape::read.tree(text = nwk)
    }
    return(tree)
}
