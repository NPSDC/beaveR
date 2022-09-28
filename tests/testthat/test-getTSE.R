test_that("check inputs", {
    skip("skip")
  expect_error(getTSE())
  expect_error(getTSE("invalid", "quant file"))

  dir <- "/home/noor/cbcb_rob/Uncertainity/brain_sim_nodtu/mode=gc_bias/post_type=gibbs_nrep=100_tf=100"
  clust_file <- file.path(dir, "terminus/no_threshold0/cluster_nwk.txt")
  expect_error(getTSE(clust_file, "ff"))

  quant_dir <- file.path(dir, "out_sal")
  samples <- as.vector(outer(c(1), c(1,2), function(x,y) paste(x,y,sep="_")))
  quant_files <- file.path(quant_dir, samples, "quant.sf")
  coldata <- data.frame(files=quant_files)
  expect_error(getTSE(clust_file, coldata))

  coldata <- data.frame(files=quant_files, names=samples, condition = as.factor(rep(c(1,2),each=1)))
  # expect_message(getTSE(clust_file, coldata))
  #
  clust_file <-"test-getTSE.R"
  expect_error(getTSE(clust_file, coldata))

  clust_file <- file.path(dir, "terminus/no_threshold0/cluster_nwk.txt")
  expect_message(getTSE(clust_file, coldata))
})

# test_that("wrong Term file", {
#
# })
