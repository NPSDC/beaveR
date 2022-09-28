test_that("rowAgg with no groups", {
    mat <- matrix(paste("a", seq(20)), nrow=4, ncol=5)
    expect_error(performRowSumAgg(mat))
    mat <- matrix(seq(20), nrow=4, ncol=5)
    aggMat <- performRowSumAgg(mat)
    expect_equal(mat, aggMat)

    expect_error(performRowSumAgg(mat[1,]))
    expect_error(performRowSumAgg(mat[,1]))

    # dir <- "/home/noor/cbcb_rob/Uncertainity/brain_sim_nodtu/mode=gc_bias/post_type=gibbs_nrep=100_tf=100"
    # quant_dir <- file.path(dir, "out_sal")
    # samples <- as.vector(outer(c(1), c(1,2), function(x,y) paste(x,y,sep="_")))
    # quant_files <- file.path(quant_dir, samples, "quant.sf")
    # coldata <- data.frame(files=quant_files, names=samples, condition = as.factor(rep(c(1,2),each=1)))

    # se <- tximeta::tximeta(coldata)
})

test_that("rowAgg with rowInds", {
    rowInds <- list(c(1:2),c(3:5))
    mat <- matrix(paste("a", seq(20)), nrow=4, ncol=5)
    expect_error(performRowSumAgg(mat))
    mat <- matrix(seq(20), nrow=4, ncol=5)

    expect_error(performRowSumAgg(mat, rowInds)) ## rowInds larger than matrix
    expect_error(performRowSumAgg(mat[1,], rowInds)) ## rowInds larger than matrix
    expect_error(performRowSumAgg(mat[1,], c(1:3))) ## rowInds larger than matrix

    rowInds <- list(c(1:2),c(3:4))
    aggMat <- performRowSumAgg(mat, rowInds)
    expected <- matrix(seq(3,39,4), nrow=2, ncol=5)
    dimnames(expected)=list(c("Group1", "Group2"),c())
    expect_equal(aggMat, expected)

    rowInds <- list("g1"=c(1:3),"g2"=c(4))
    colnames(mat) <- paste("a", c(1:5), sep="")
    aggMat <- performRowSumAgg(mat, rowInds)
    expected <- matrix(c(6,4,18,8,30,12,42,16,54,20), nrow=2, ncol=5)
    dimnames(expected)=list(c("g1", "g2"),colnames(mat))
    expect_equal(aggMat, expected)

    ### Sparse Matrix
    aggMat <- performRowSumAgg(Matrix::Matrix(mat,sparse=T), rowInds)
    expected <- Matrix::Matrix(c(6,4,18,8,30,12,42,16,54,20), nrow=2, ncol=5, sparse=T)
    dimnames(expected)=list(c("g1", "g2"),colnames(mat))
    expect_equal(aggMat, expected)
    #
    #
    # expect_equal(mat[1,], aggMat)
    #
    # aggMat <- performRowSumAgg(mat[,1])
    # expect_equal(mat[,1], aggMat)
    #
    # aggMat <- performRowSumAgg(mat[,1])
    # dir <- "/home/noor/cbcb_rob/Uncertainity/brain_sim_nodtu/mode=gc_bias/post_type=gibbs_nrep=100_tf=100"
    # quant_dir <- file.path(dir, "out_sal")
    # samples <- as.vector(outer(c(1), c(1,2), function(x,y) paste(x,y,sep="_")))
    # quant_files <- file.path(quant_dir, samples, "quant.sf")
    # coldata <- data.frame(files=quant_files, names=samples, condition = as.factor(rep(c(1,2),each=1)))

    # se <- tximeta::tximeta(coldata)
})
