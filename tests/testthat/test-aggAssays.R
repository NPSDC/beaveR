test_that("aggAssays", {
    dir <- "/home/noor/cbcb_rob/Uncertainity/brain_sim_nodtu/mode=gc_bias/post_type=gibbs_nrep=100_tf=100"
    clustFile <- file.path(dir, "terminus/no_threshold0/cluster_nwk.txt")
    quantDir <- file.path(dir, "out_sal")
    samples <- as.vector(outer(c(1), c(1,2), function(x,y) paste(x,y,sep="_")))
    quantFiles <- file.path(quantDir, samples, "quant.sf")
    coldata <- data.frame(files=quantFiles, names=samples, condition = as.factor(rep(c(1,2),each=1)))
    trees <- ape::read.tree(clustFile)
    expect_error(aggAssays(trees, trees))

    # se <- tximeta::tximeta(coldata)
    expect_error(aggAssays(trees[[1]], trees))

    # tree <- mergeTrees(trees)
    # tree <- mergeTrees(trees, tnames=rownames(se))

})
