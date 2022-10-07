test_that("findNodeInformation works", {
    dir <- system.file("extdata", package="beaveR")
    dir <- file.path(dir, "brain_sim_nodtu_small_example")
    samples <- as.vector(outer(c(1:6), c(1,2), function(x,y) paste(x,y,sep="_")))
    quantDir <- file.path(dir, "out_sal")
    quantFiles <- file.path(quantDir, samples, "quant.sf")
    clustFile <- file.path(dir, "cluster_nwk.txt")
    coldata <- data.frame(files = quantFiles, names = samples)
    tse <- buildTSE(clustFile, coldata)
    expect_error(findNodeInformation(coldata, coldata))
    expect_error(findNodeInformation(tse, "aa"))
    expect_error(findNodeInformation(tse, 400))
    expect_error(findNodeInformation(tse, 10, "aa"))
    df <- findNodeInformation(tse, 100)
    expect_equal(nrow(df), 1)
    df <- findNodeInformation(tse, 100, "children")
    expect_equal(nrow(df), 1)

    df <- findNodeInformation(tse, 300, "tips")
    expect_equal(nrow(df), 7)

    df <- findNodeInformation(tse, 300, "all")
    expect_equal(nrow(df), 11)

    df <- findNodeInformation(tse, 300)
    expect_equal(nrow(df), 1)

    expect_equal(ncol(df), 4)
})
