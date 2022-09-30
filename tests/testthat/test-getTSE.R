test_that("check inputs", {
  expect_error(getTSE())
  expect_error(getTSE("invalid", "quant file"))

  dir <- "../../extdata/brain_sim_nodtu_small_example"
  clustFile <- file.path(dir, "group_nwk.txt")
  expect_error(getTSE(clustFile, "ff"))

  quantDir <- file.path(dir, "out_sal")
  samples <- as.vector(outer(seq(6), c(1,2), function(x,y) paste(x,y,sep="_")))
  quantFiles <- file.path(quantDir, samples, "quant.sf")
  coldata <- data.frame(files=quantFiles)
  expect_error(getTSE(clustFile, coldata))

  coldata <- data.frame(files=quantFiles, names=samples)
  # expect_message(getTSE(clust_file, coldata))
  #
  clustFile <-"test-getTSE.R"
  expect_error(getTSE(clustFile, coldata))

  clustFile <- file.path(dir, "group_nwk.txt")
  tse <- getTSE(clustFile, coldata)
  expect_s4_class(tse, "TreeSummarizedExperiment")
})

# test_that("wrong Term file", {
#
# })
