test_that("Test trainRankModel", {
    seu_sc <- readRDS(system.file("extdata", "seu_sc.rds",
        package = "RankMap"
    ))

    mat <- extractData(seu_sc)

    model <- trainRankModel(mat, seu_sc$cell_type)
    expect_true(inherits(model, "glmnet"))

    model <- trainRankModel(
        data = mat,
        labels = seu_sc$cell_type,
        cv = TRUE,
        nfolds = 3
    )
    expect_true(inherits(model, "cv.glmnet"))
})


test_that("Test predictRankModel", {
    seu_sc <- readRDS(system.file("extdata", "seu_sc.rds",
        package = "RankMap"
    ))

    seu_xen <- readRDS(system.file("extdata", "seu_xen.rds",
        package = "RankMap"
    ))

    common_genes <- intersect(rownames(seu_sc), rownames(seu_xen))
    mat <- extractData(seu_sc)[common_genes, ]
    new_mat <- extractData(seu_xen)[common_genes, ]

    model <- trainRankModel(mat, seu_sc$cell_type)

    pred <- predictRankModel(model, new_mat)

    expect_true(is.character(pred))

    pred <- predictRankModel(model, new_mat, return_confidence = TRUE)

    expect_true(is.data.frame(pred))

    pred <- predictRankModel(model, new_mat, return_probs = TRUE)

    expect_true(is.matrix(pred))

    expect_error(predictRankModel(model, new_mat,
        return_confidence = TRUE,
        return_probs = TRUE
    ))

    model <- trainRankModel(
        data = mat,
        labels = seu_sc$cell_type,
        cv = TRUE,
        nfolds = 3
    )
    pred <- predictRankModel(model, new_mat)
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
