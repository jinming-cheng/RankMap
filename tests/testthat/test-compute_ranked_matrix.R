test_that("Test functions in compute_ranked_matrix", {
    seu_sc <- readRDS(system.file("extdata", "seu_sc.rds",
        package = "RankMap"
    ))

    mat <- ExtractData(seu_sc)

    expect_true(is.matrix(as.matrix(mat)))

    expect_error(ExtractData(as.data.frame(mat)))

    seu_sc <- Seurat::ScaleData(seu_sc)
    sce <- Seurat::as.SingleCellExperiment(seu_sc)

    mat <- ExtractData(sce)

    expect_true(is.matrix(as.matrix(mat)))

    masked <- MaskTopKGenes(mat)

    expect_true(is.matrix(masked))

    mat_ranked <- ComputeRankedMatrix(mat)

    expect_true(is.matrix(mat_ranked))

    mat_ranked <- ComputeRankedMatrix(mat, use_data = TRUE)

    expect_equal(mat, mat_ranked)
})
