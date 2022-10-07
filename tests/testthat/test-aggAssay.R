test_that("rowAgg without rowInds", {
    mat <- matrix(paste("a", seq(20)), nrow = 4, ncol = 5)
    expect_error(performRowSumAgg(mat))
    mat <- matrix(seq(20), nrow = 4, ncol = 5)
    aggMat <- performRowSumAgg(mat)
    expect_equal(mat, aggMat)

    expect_error(performRowSumAgg(mat[1, ]))
    expect_error(performRowSumAgg(mat[, 1]))

})

test_that("rowAgg with rowInds", {

    rowInds <- list(c(1:2), c(3:5))
    mat <- matrix(paste("a", seq(20)), nrow = 4, ncol = 5)
    expect_error(performRowSumAgg(mat))
    mat <- matrix(seq(20), nrow = 4, ncol = 5)

    ## rowInds larger than rows of the matrix
    expect_error(performRowSumAgg(mat, rowInds))
    expect_error(performRowSumAgg(mat[1, ], rowInds))

    ## rowInds not list
    expect_error(performRowSumAgg(mat[1, ], c(1:3)))

    rowInds <- list(c(1:2), c(3:4))
    aggMat <- performRowSumAgg(mat, rowInds)
    expected <- matrix(seq(3, 39, 4), nrow = 2, ncol = 5)
    dimnames(expected) = list(c("Group1", "Group2"), c())
    expect_equal(aggMat, expected)

    rowInds <- list(g1 = c(1:3), g2 = c(4))
    colnames(mat) <- paste("a", c(1:5), sep = "")
    aggMat <- performRowSumAgg(mat, rowInds)
    expected <- matrix(c(6, 4, 18, 8, 30, 12, 42, 16, 54, 20), nrow = 2, ncol = 5)
    dimnames(expected) = list(c("g1", "g2"), colnames(mat))
    expect_equal(aggMat, expected)

    ### Sparse Matrix
    aggMat <- performRowSumAgg(Matrix::Matrix(mat, sparse = T), rowInds)
    expected <- Matrix::Matrix(c(6, 4, 18, 8, 30, 12, 42, 16, 54, 20), nrow = 2, ncol = 5,
        sparse = T)
    dimnames(expected) = list(c("g1", "g2"), colnames(mat))
    expect_equal(aggMat, expected)
})

test_that("colAgg without colInds", {
    mat <- matrix(paste("a", seq(20)), nrow = 4, ncol = 5)
    expect_error(performRowSumAgg(mat))
    mat <- matrix(seq(20), nrow = 4, ncol = 5)
    aggMat <- performColMeanAgg(mat)
    expect_equal(mat, aggMat)

    expect_error(performColMeanAgg(mat[1, ]))
    expect_error(performColMeanAgg(mat[, 1]))

})

test_that("colAgg with colInds", {
    colInds <- list(c(1:2), c(3:6))
    mat <- matrix(paste("a", seq(20)), nrow = 4, ncol = 5)
    expect_error(performColMeanAgg(mat))
    mat <- matrix(seq(20), nrow = 4, ncol = 5)

    ## colInds larger than number of columns
    expect_error(performColMeanAgg(mat, colInds))
    expect_error(performColMeanAgg(mat[1, ], rowInds))

    ## colInds not list
    expect_error(performColMeanAgg(mat[1, ], c(1:3)))

    colInds <- list(c(1:2), c(3:5))
    aggMat <- performColMeanAgg(mat, colInds)
    expected <- matrix(c(3, 4, 5, 6, 13, 14, 15, 16), nrow = 4, ncol = 2)
    dimnames(expected) = list(c(), c("Group1", "Group2"))
    expect_equal(aggMat, expected)

    colInds <- list(g1 = c(1:4), g2 = c(5))
    rownames(mat) <- paste("a", c(1:4), sep = "")
    aggMat <- performColMeanAgg(mat, colInds)
    expected <- matrix(c(7, 8, 9, 10, 17, 18, 19, 20), nrow = 4, ncol = 2)
    dimnames(expected) = list(rownames(mat), c("g1", "g2"))
    expect_equal(aggMat, expected)

    ### Sparse Matrix
    aggMat <- performColMeanAgg(Matrix::Matrix(mat, sparse = T), colInds)
    expected <- Matrix::Matrix(c(7, 8, 9, 10, 17, 18, 19, 20), nrow = 4, ncol = 2, sparse = T)
    dimnames(expected) = list(rownames(mat), c("g1", "g2"))
    expect_equal(aggMat, expected)

})

test_that("aggAssay", {
    mat <- matrix(seq(20), nrow = 5, ncol = 4)
    tree <- mat
    expect_error(aggAssay(tree, 5, mat))
    tree <- ape::rtree(5)
    expect_error(aggAssay(tree, "aa", mat))
    expect_error(aggAssay(tree, 5, mat[1, ]))
    expect_error(aggAssay(tree, 5, mat[1, ]))
    expect_error(aggAssay(tree, c(1:4, 10), mat))
    expect_error(aggAssay(tree, c(1:9), mat, list(1:5)))
    aMat <- aggAssay(tree, c(1:5), mat)
    expect_equal(aMat, mat)
    aMat <- aggAssay(tree, c(1:6), mat)
    eMat <- rbind(mat, colSums(mat))
    rownames(eMat) <- c(rep("", 5), "Node6")
    expect_equal(aMat, eMat)
    aMat <- aggAssay(tree, c(6), mat)
    eMat <- matrix(colSums(mat), nrow = 1, ncol = 4)
    rownames(eMat) <- "Node6"
    expect_equal(aMat, eMat)
    aMat <- aggAssay(tree, c(2), mat)
    expect_equal(aMat, matrix(mat[2, ], nrow = 1))
})
