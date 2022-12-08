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

    l <- length(tree$tip)
    expect_s4_class(seAgg, "SummarizedExperiment")
    expect_equal(SummarizedExperiment::assays(seAgg)[["counts"]][l + 1, ], colSums(SummarizedExperiment::assays(seAgg)[["counts"]][1:l,
        ]))
    expect_equal(SummarizedExperiment::assays(seAgg)[["infRep10"]][l + 1, ], colSums(SummarizedExperiment::assays(seAgg)[["infRep10"]][1:l,
        ]))
    expect_equal(SummarizedExperiment::assays(seAgg)[["infRep1"]][l + 1, ], colSums(SummarizedExperiment::assays(seAgg)[["infRep1"]][1:l,
        ]))
    expect_equal(SummarizedExperiment::assays(seAgg)[["counts"]][l + 2, ], colSums(SummarizedExperiment::assays(seAgg)[["counts"]][phangorn::Descendants(tree,
        l + 2)[[1]], ]))
    expect_equal(SummarizedExperiment::assays(seAgg)[["counts"]][l + 2, ], colSums(SummarizedExperiment::assays(seAgg)[["counts"]][phangorn::Descendants(tree,
                                                                                                                                                         l + 2)[[1]], ]))
    expect_equal(nrow(seAgg), length(tree$tip) + tree$Nnode)
    expect_equal(SummarizedExperiment::assays(seAgg)[["length"]][l - 1, ], SummarizedExperiment::assays(seAgg)[["length"]][l - 1,])
})
