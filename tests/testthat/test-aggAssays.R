test_that("aggAssays", {
  dir <- system.file("extdata", package = "beaveR")
  dir <- file.path(dir, "brain_sim_nodtu_small_example")
  clustFile <- file.path(dir, "cluster_nwk.txt")
  quantDir <- file.path(dir, "out_sal")
  samples <- as.vector(outer(c(1:6), c(1, 2), function(x, y) paste(x, y, sep = "_")))
  quantFiles <- file.path(quantDir, samples, "quant.sf")
  coldata <- data.frame(files = quantFiles, names = samples)
  trees <- ape::read.tree(clustFile)
  expect_error(aggAssays(trees, trees))

  se <- tximeta::tximeta(coldata)
  tree <- mergeTrees(trees, tnames = rownames(se))
  txpsFilt <- getTxps(txps = NULL, y = se, minN = 1)
  tree <- mergeTreeWithSE(tree, txpsFilt)
  seAgg <- aggAssays(tree, se[tree$tip, ])

  ## Checking for length
  l <- length(tree$tip)
  expect_equal(assays(se)[["length"]], assays(seAgg)[["length"]][1:l,])
  ## Tree nodes that completely contain genes -
  ## checking that our length procedure is exact as tximeta
  gs <- tximeta::summarizeToGene(se)
  gInds <- lapply(SummarizedExperiment::rowData(gs)[["tx_ids"]], function(g) match(unlist(g), tree$tip.label)) ## index of transcripts belonging to genes in the tree
  gInds <- gInds[sapply(gInds, function(x) sum(is.na(x)) == 0)] ## Only those genes whose entire transcript matches the tree transcripts
  # ENSG00000221843.4 - 374
  # ENSG00000240524.2 - 6
  # ENSG00000188566.13 - 375
  # ENSG00000156395.12 - 227
  # ENSG00000143344.15 - 290
  expect_equal(assays(seAgg)[["length"]][290,], assays(gs)[["length"]]["ENSG00000143344.15",])
  expect_equal(assays(seAgg)[["length"]][227,], assays(gs)[["length"]]["ENSG00000156395.12",])
  expect_equal(assays(seAgg)[["length"]][375,], assays(gs)[["length"]]["ENSG00000188566.13",])
  expect_equal(assays(seAgg)[["length"]][6,], assays(gs)[["length"]]["ENSG00000240524.2",])
  expect_equal(assays(seAgg)[["length"]][374,], assays(gs)[["length"]]["ENSG00000221843.4",])

  expect_s4_class(seAgg, "SummarizedExperiment")
  expect_equal(SummarizedExperiment::assays(seAgg)[["counts"]][l + 1, ], colSums(SummarizedExperiment::assays(seAgg)[["counts"]][1:l, ]))
  expect_equal(SummarizedExperiment::assays(seAgg)[["infRep10"]][l + 1, ], colSums(SummarizedExperiment::assays(seAgg)[["infRep10"]][1:l, ]))
  expect_equal(SummarizedExperiment::assays(seAgg)[["infRep1"]][l + 1, ], colSums(SummarizedExperiment::assays(seAgg)[["infRep1"]][1:l, ]))
  expect_equal(SummarizedExperiment::assays(seAgg)[["counts"]][l + 2, ], colSums(SummarizedExperiment::assays(seAgg)[["counts"]][phangorn::Descendants(
    tree,
    l + 2
  )[[1]], ]))
  expect_equal(SummarizedExperiment::assays(seAgg)[["counts"]][l + 2, ], colSums(SummarizedExperiment::assays(seAgg)[["counts"]][phangorn::Descendants(
    tree,
    l + 2
  )[[1]], ]))
  expect_equal(nrow(seAgg), length(tree$tip) + tree$Nnode)
  expect_equal(SummarizedExperiment::assays(seAgg)[["length"]][l - 1, ], SummarizedExperiment::assays(seAgg)[["length"]][l - 1, ])
})
