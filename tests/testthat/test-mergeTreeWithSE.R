test_that("mergeTreeWithSE working", {
  dir <- system.file("extdata", package = "beaveR")
  samples <- as.vector(outer(c(1:6), c(1, 2), function(x, y) paste(x, y, sep = "_")))
  clustFile <- file.path(dir, "brain_sim_nodtu_small_example", "cluster_nwk.txt")
  quantFiles <- file.path(dir, "brain_sim_nodtu_small_example", "out_sal", samples, "quant.sf")
  trees <- ape::read.tree(clustFile)

  # quantDir <- file.path(dir, )


  coldata <- data.frame(files = quantFiles, names = samples)
  se <- tximeta::tximeta(coldata)

  tree <- mergeTrees(trees)
  expect_error(mergeTreeWithSE(tree, se))

  tree <- mergeTrees(trees, tnames = rownames(se))
  txpsFilt <- getTxps(txps = NULL, y = se, minN = 1)
  tree <- mergeTreeWithSE(tree, txpsFilt)
  expect_true(all(tree$tip %in% rownames(se)))
})
