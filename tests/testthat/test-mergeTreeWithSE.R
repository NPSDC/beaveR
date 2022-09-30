test_that("mergeTreeWithSE working", {
    dir <- "../../extdata/brain_sim_nodtu_small_example"
    clustFile <- file.path(dir, "group_nwk.txt")
    trees <- ape::read.tree(clustFile)

    quantDir <- file.path(dir, "out_sal")
    samples <- as.vector(outer(c(1:6), c(1,2), function(x,y) paste(x,y,sep="_")))
    quantFiles <- file.path(quantDir, samples, "quant.sf")
    coldata <- data.frame(files=quantFiles, names=samples)
    se <- tximeta::tximeta(coldata)

    tree <- mergeTrees(trees)
    expect_error(mergeTreeWithSE(tree, se))

    tree <- mergeTrees(trees, tnames=rownames(se))
    txpsFilt <- getTxps(txps=NULL, y=se, minN=1)
    tree <- mergeTreeWithSE(tree, txpsFilt)
    expect_true(all(tree$tip %in% rownames(se)))
})
