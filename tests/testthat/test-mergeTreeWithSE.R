test_that("mergeTreeWithSE working", {
    skip("skip")
    dir <- "/home/noor/cbcb_rob/Uncertainity/brain_sim_nodtu/mode=gc_bias/post_type=gibbs_nrep=100_tf=100"
    clust_file <- file.path(dir, "terminus/no_threshold0/cluster_nwk.txt")
    trees <- ape::read.tree(clust_file)

    quant_dir <- file.path(dir, "out_sal")
    samples <- as.vector(outer(c(1), c(1,2), function(x,y) paste(x,y,sep="_")))
    quant_files <- file.path(quant_dir, samples, "quant.sf")
    coldata <- data.frame(files=quant_files, names=samples, condition = as.factor(rep(c(1,2),each=1)))
    se <- tximeta::tximeta(coldata)

    tree <- mergeTrees(trees)
    expect_error(mergeTreeWithSE(tree, se))

    tree <- mergeTrees(trees, tnames=rownames(se))
    txps <- getTxps(txps=NULL, y=se, minN=1)
    tree <- mergeTreeWithSE(tree, se[txps,])
    expect_true(all(tree$tip %in% rownames(se)))
})
