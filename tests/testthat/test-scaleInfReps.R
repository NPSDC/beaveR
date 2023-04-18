# test_that("compSizeFactors", {
#     dir <- system.file("extdata", package = "beaveR")
#     dir <- file.path(dir, "brain_sim_nodtu_small_example")
#     clustFile <- file.path(dir, "cluster_nwk.txt")
#     expect_error(computeSizeFactors(clustFile)) ##TSE
#
#     quantDir <- file.path(dir, "out_sal")
#     samples <- as.vector(outer(seq(6), c(1, 2), function(x, y) paste(x, y, sep = "_")))
#     quantFiles <- file.path(quantDir, samples, "quant.sf")
#     coldata <- data.frame(files = quantFiles, names = samples)
#     tse <- buildTSE(clustFile, coldata)
#     expect_error(computeSizeFactors(clustFile)) ##TSE
#
#     expect_error(computeSizeFactors(tse, type="dd")) ##type
#
#     tse <- computeSizeFactors(tse)
#     sizeMat <- metadata(tse)[["sf"]]
#
#     expect_identical(dim(sizeMat), as.integer(c(100,12)))
#     expect_equal(metadata(tse)[["type"]], "txp")
#     expect_error(computeSizeFactors(tse))
#     tse <- computeSizeFactors(tse, force = TRUE)
#     expect_equal(metadata(tse)[["type"]], "txp")
#     tse <- buildTSE(clustFile, coldata)
#     tse <- computeSizeFactors(tse, type = "gene")
#     expect_equal(metadata(tse)[["type"]], "gene")
# })

test_that("scaleInfReps", {
  dir <- system.file("extdata", package = "beaveR")
  dir <- file.path(dir, "brain_sim_nodtu_small_example")
  clustFile <- file.path(dir, "cluster_nwk.txt")
  expect_error(scaleInfReps(clustFile)) ## TSE

  quantDir <- file.path(dir, "out_sal")
  samples <- as.vector(outer(seq(6), c(1, 2), function(x, y) paste(x, y, sep = "_")))
  quantFiles <- file.path(quantDir, samples, "quant.sf")
  coldata <- data.frame(files = quantFiles, names = samples)
  tse <- buildTSE(clustFile, coldata)
  tseDup <- tse

  expect_error(scaleInfReps(tse))
  expect_error(scaleInfReps(tse, "aaa")) ## sf

  # print(tail(assays(tse)[["infRep2"]]))
  sf <- matrix(1, nrow = 10, ncol = 12)
  expect_error(scaleInfReps(tse, szMat = sf))
  sf <- matrix(1, nrow = 100, ncol = 12)
  tse <- scaleInfReps(tse, szMat = sf)
  # print(tail(assays(tse)[["infRep2"]]))

  expect_true(all(sf == metadata(tse)[["sf"]]))
  expect_error(scaleInfReps(tse))
  expect_error(scaleInfReps(tse, sf))

  tse <- computeSizeFactors(tse, force = TRUE)
  expect_error(scaleInfReps(tse))
  tse <- scaleInfReps(tse, force = TRUE)

  expect_error(scaleInfReps(tse, sf))
  expect_identical(dim(tse), as.integer(c(387, 12)))

  tse <- buildTSE(clustFile, coldata)
  tse <- computeSizeFactors(tse)

  tse <- scaleInfReps(tse)
  expect_identical(dim(tse), as.integer(c(387, 12)))

  tseDup <- fishpond::scaleInfReps(tseDup[1:203, ])
  expect_equal(assays(tse)[["infRep10"]][1:203, ], assays(tseDup)[["infRep10"]])
})
