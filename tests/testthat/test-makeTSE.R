test_that("check inputs", {
    expect_error(makeTSE())
    expect_error(makeTSE("invalid", "quant file"))

    dir <- "../../extdata/brain_sim_nodtu_small_example"
    clustFile <- file.path(dir, "group_nwk.txt")
    expect_error(makeTSE(clustFile, "ff"))

    quantDir <- file.path(dir, "out_sal")
    samples <-
        as.vector(outer(seq(6), c(1, 2), function(x, y)
            paste(x, y, sep = "_")))
    quantFiles <- file.path(quantDir, samples, "quant.sf")
    coldata <- data.frame(files = quantFiles)
    expect_error(makeTSE(clustFile, coldata))

    coldata <- data.frame(files = quantFiles, names = samples)
    # expect_message(makeTSE(clust_file, coldata))
    #
    clustFile <- "test-makeTSE.R"
    expect_error(makeTSE(clustFile, coldata))

    clustFile <- file.path(dir, "group_nwk.txt")
    tse <- makeTSE(clustFile, coldata)
    expect_s4_class(tse, "TreeSummarizedExperiment")
    mirv <- SummarizedExperiment::mcols(tse)[["meanInfRV"]]
    expect_equal(length(mirv), nrow(tse))
})

# test_that("wrong Term file", {
#
# })
