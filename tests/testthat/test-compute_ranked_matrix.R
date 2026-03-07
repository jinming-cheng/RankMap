test_that("Test functions in compute_ranked_matrix", {
    seu_sc <- readRDS(system.file("extdata", "seu_sc.rds",
        package = "RankMap"
    ))

    mat <- extractData(seu_sc)

    expect_true(is.matrix(as.matrix(mat)))

    expect_error(extractData(as.data.frame(mat)))

    seu_sc <- Seurat::ScaleData(seu_sc)
    sce <- Seurat::as.SingleCellExperiment(seu_sc)

    mat <- extractData(sce)

    expect_true(is.matrix(as.matrix(mat)))

    masked <- maskTopKGenes(mat)

    expect_true(is.matrix(masked))

    mat_ranked <- computeRankedMatrix(mat)

    expect_true(is.matrix(mat_ranked))

    mat_ranked <- computeRankedMatrix(mat, use_data = TRUE)

    expect_equal(mat, mat_ranked)
})
