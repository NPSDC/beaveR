test_that("testSolveObj", {
    dir <- "../../extdata/brain_sim_nodtu_small_example"
    clustFile <- file.path(dir, "group_nwk.txt")
    quantDir <- file.path(dir, "out_sal")
    samples <- as.vector(outer(seq(6), c(1,2), function(x,y) paste(x,y,sep="_")))
    quantFiles <- file.path(quantDir, samples, "quant.sf")
    coldata <- data.frame(files=quantFiles, names=samples)
    tse <- getTSE(clustFile, coldata)
    tree <- TreeSummarizedExperiment::rowTree(tse)
    gamma <- 0.1
    l <- length(tree$tip)
    sizeDesc <- sapply(phangorn::Descendants(tree, c(1:l,l+1:tree$Nnode)),length)
    metVec <- SummarizedExperiment::mcols(tse)[["meanInfRV"]] + gamma*ape::node.depth(tree, 2)
    metVec <- metVec*sizeDesc
    expect_error(solveObj(coldata, coldata))
    expect_error(solveObj(tse, metVec[1:10]))
    expect_error(solveObj(tse, metVec, "ss"))

    opt <- solveObj(tse, metVec, "min")
    print(length(opt[["cut"]]))
    print(length(opt[["optVal"]]))
})
