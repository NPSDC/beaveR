test_that("transcript NULL", {
  dir <- system.file("extdata", package = "beaveR")
  dir <- file.path(dir, "brain_sim_nodtu_small_example")
  clustFile <- file.path(dir, "cluster_nwk.txt")
  quantDir <- file.path(dir, "out_sal")
  samples <- as.vector(outer(seq(6), c(1, 2), function(x, y) paste(x, y, sep = "_")))
  quantFiles <- file.path(quantDir, samples, "quant.sf")
  coldata <- data.frame(files = quantFiles, names = samples, condition = as.factor(rep(c(
    1,
    2
  ), each = 1)))
  se <- tximeta::tximeta(coldata)
  expect_error(getTxps(rownames(se)))
  expect_error(getTxps(1, rownames(se)))
  txps <- getTxps(rownames(se), se)
  expect_true(all(txps %in% rownames(se)))
  expect_error(getTxps(txps = NULL, se, dd = 1e+10))
  txps <- getTxps(txps = NULL, se, minN = 1)
  expect_true(all(txps %in% rownames(se)))
})
