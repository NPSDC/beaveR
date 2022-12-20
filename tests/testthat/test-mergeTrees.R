test_that("mergeTree working", {
  dir <- system.file("extdata", package = "beaveR")
  dir <- file.path(dir, "brain_sim_nodtu_small_example")
  clustFile <- file.path(dir, "cluster_nwk.txt")
  trees <- ape::read.tree(clustFile)

  tree <- mergeTrees(trees)
  lengths <- sum(sapply(trees, function(tree) length(tree$tip)))
  expect_equal(lengths, length(tree$tip))

  quantDir <- file.path(dir, "out_sal")
  samples <- as.vector(outer(c(1:6), c(1, 2), function(x, y) paste(x, y, sep = "_")))
  quantFiles <- file.path(quantDir, samples, "quant.sf")
  coldata <- data.frame(files = quantFiles, names = samples)
  se <- tximeta::tximeta(coldata)

  expect_error(mergeTrees(trees, tnames = rownames(se)[1:100]))
  expect_error(mergeTrees(trees, tnames = 1:10))

  treesCopy <- trees
  treesCopy[[1]]$tip.label[1] <- as.character(as.numeric(treesCopy[[1]]$tip.label[1]) +
    1e+05)
  expect_error(mergeTrees(treesCopy, tnames = rownames(se)))

  tree <- mergeTrees(trees, tnames = rownames(se))
  expect_equal(lengths, length(tree$tip))
})
