test_that("transcript NULL", {
    skip('skip')
    dir <- "/home/noor/cbcb_rob/Uncertainity/brain_sim_nodtu/mode=gc_bias/post_type=gibbs_nrep=100_tf=100"
    clust_file <- file.path(dir, "terminus/no_threshold0/cluster_nwk.txt")
    quant_dir <- file.path(dir, "out_sal")
    samples <- as.vector(outer(c(1), c(1,2), function(x,y) paste(x,y,sep="_")))
    quant_files <- file.path(quant_dir, samples, "quant.sf")
    coldata <- data.frame(files=quant_files, names=samples, condition = as.factor(rep(c(1,2),each=1)))
    se <- tximeta::tximeta(coldata)
    expect_error(getTxps(rownames(se)))
    expect_error(getTxps(1,rownames(se)))
    txps <- getTxps(rownames(se), se)
    expect_true(all(txps %in% rownames(se)))
    expect_error(getTxps(txps=NULL, se))
    expect_error(getTxps(txps=NULL, se, minCount=1e10))
    expect_error(getTxps(txps=NULL, se, dd=1e10))
    txps <- getTxps(txps=NULL, se, minN=1)
    expect_true(all(txps %in% rownames(se)))
})