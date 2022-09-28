test_that("mergeTree working", {
    skip("skip")
    dir <- "/home/noor/cbcb_rob/Uncertainity/brain_sim_nodtu/mode=gc_bias/post_type=gibbs_nrep=100_tf=100"
    clust_file <- file.path(dir, "terminus/no_threshold0/cluster_nwk.txt")
    trees <- ape::read.tree(clust_file)

    tree <- mergeTrees(trees)
    lengths <- sum(sapply(trees, function(tree) length(tree$tip)))
    expect_equal(lengths, length(tree$tip))

    quant_dir <- file.path(dir, "out_sal")
    samples <- as.vector(outer(c(1), c(1,2), function(x,y) paste(x,y,sep="_")))
    quant_files <- file.path(quant_dir, samples, "quant.sf")
    coldata <- data.frame(files=quant_files, names=samples, condition = as.factor(rep(c(1,2),each=1)))
    se <- tximeta::tximeta(coldata)
    expect_error(mergeTrees(trees, tnames = rownames(se)[1:100]))
    expect_error(mergeTrees(trees, tnames = 1:10))
    tree <- mergeTrees(trees, tnames = rownames(se))
    expect_equal(lengths, length(tree$tip))
})
