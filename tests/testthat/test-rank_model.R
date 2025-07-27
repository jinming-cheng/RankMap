test_that("Test TrainRankModel", {
    seu_sc <- readRDS(system.file("extdata", "seu_sc.rds",
        package = "RankMap"
    ))

    mat <- ExtractData(seu_sc)

    model <- TrainRankModel(mat, seu_sc$cell_type)
    expect_true(inherits(model, "glmnet"))

    model <- TrainRankModel(
        data = mat,
        labels = seu_sc$cell_type,
        cv = TRUE,
        nfolds = 3
    )
    expect_true(inherits(model, "cv.glmnet"))
})


test_that("Test PredictRankModel", {
    seu_sc <- readRDS(system.file("extdata", "seu_sc.rds",
        package = "RankMap"
    ))

    seu_xen <- readRDS(system.file("extdata", "seu_xen.rds",
        package = "RankMap"
    ))

    common_genes <- intersect(rownames(seu_sc), rownames(seu_xen))
    mat <- ExtractData(seu_sc)[common_genes, ]
    new_mat <- ExtractData(seu_xen)[common_genes, ]

    model <- TrainRankModel(mat, seu_sc$cell_type)

    pred <- PredictRankModel(model, new_mat)

    expect_true(is.character(pred))

    pred <- PredictRankModel(model, new_mat, return_confidence = TRUE)

    expect_true(is.data.frame(pred))

    pred <- PredictRankModel(model, new_mat, return_probs = TRUE)

    expect_true(is.matrix(pred))

    expect_error(PredictRankModel(model, new_mat,
        return_confidence = TRUE,
        return_probs = TRUE
    ))

    model <- TrainRankModel(
        data = mat,
        labels = seu_sc$cell_type,
        cv = TRUE,
        nfolds = 3
    )
    pred <- PredictRankModel(model, new_mat)
    expect_true(is.character(pred))
})


test_that("Test RankMap", {
    seu_sc <- readRDS(system.file("extdata", "seu_sc.rds",
        package = "RankMap"
    ))

    seu_xen <- readRDS(system.file("extdata", "seu_xen.rds",
        package = "RankMap"
    ))

    pred_df <- RankMap(
        ref_data = seu_sc,
        ref_labels = seu_sc$cell_type,
        new_data = seu_xen,
        threshold = 0.5
    )
    expect_true(is.data.frame(pred_df))

    pred_list <- RankMap(
        ref_data = seu_sc,
        ref_labels = seu_sc$cell_type,
        new_data = seu_xen,
        return_model = TRUE
    )
    expect_true(is.list(pred_list))

    expect_error(RankMap(
        ref_data = seu_sc,
        ref_labels = seu_sc$cell_type,
        new_data = seu_xen[1:20, ]
    ))

    expect_message(RankMap(
        ref_data = seu_sc,
        ref_labels = seu_sc$cell_type,
        new_data = seu_sc
    ))
})
