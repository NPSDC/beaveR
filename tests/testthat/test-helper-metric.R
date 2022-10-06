test_that("scaledLFC", {
    dir <- "../../extdata/brain_sim_nodtu_small_example"
    clustFile <- file.path(dir, "group_nwk.txt")
    quantDir <- file.path(dir, "out_sal")
    samples <- as.vector(outer(c(1:6), c(1,2), function(x,y) paste(x,y,sep="_")))
    quantFiles <- file.path(quantDir, samples, "quant.sf")
    coldata <- data.frame(files=quantFiles, names=samples)
    tse <- makeTSE(clustFile, coldata)
    expect_error(getScaledLFC(coldata, coldata))
    expect_error(getScaledLFC(tse, "condition"))

    coldata <- data.frame(files=quantFiles, names=samples,condition=factor(rep(1:2, each=6)))
    tse <- makeTSE(clustFile, coldata)

    lfc <- getScaledLFC(tse, "condition")
    tree <- TreeSummarizedExperiment::rowTree(tse)
    l <- length(tree$tip)
    expect_equal(length(lfc), length(tree$tip)+tree$Nnode)
    expect_equal(lfc[l+1], 0)

})