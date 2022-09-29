test_that("aggAssays", {
    dir <- "/home/noor/cbcb_rob/Uncertainity/brain_sim_nodtu/mode=gc_bias/post_type=gibbs_nrep=100_tf=100"
    clustFile <- file.path(dir, "terminus/no_threshold0/cluster_nwk.txt")
    quantDir <- file.path(dir, "out_sal")
    samples <- as.vector(outer(c(1), c(1,2), function(x,y) paste(x,y,sep="_")))
    quantFiles <- file.path(quantDir, samples, "quant.sf")
    coldata <- data.frame(files=quantFiles, names=samples, condition = as.factor(rep(c(1,2),each=1)))
    trees <- ape::read.tree(clustFile)
    expect_error(aggAssays(trees, trees))

    se <- tximeta::tximeta(coldata)
    SummarizedExperiment::assays(se) <- SummarizedExperiment::assays(se)[1:13]
    tree <- mergeTrees(trees, tnames=rownames(se))
    txps <- getTxps(txps=NULL, y=se, minN=1)
    tree <- mergeTreeWithSE(tree, se[txps,])
    seAgg <- aggAssays(tree, se[tree$tip,])

    l <- length(tree$tip)
    expect_equal(SummarizedExperiment::assays(seAgg)[["counts"]][l+1,],
                 colSums(SummarizedExperiment::assays(seAgg)[["counts"]][1:l,]))
})
