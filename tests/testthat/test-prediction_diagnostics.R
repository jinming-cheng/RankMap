test_that("Test evaluatePredictionPerformance", {
    seu_sc <- readRDS(system.file("extdata", "seu_sc.rds",
        package = "RankMap"
    ))

    seu_xen <- readRDS(system.file("extdata", "seu_xen.rds",
        package = "RankMap"
    ))

    prediction_df <- RankMap(
        ref_data = seu_sc,
        ref_labels = seu_sc$cell_type,
        new_data = seu_xen
    )

    performance <- evaluatePredictionPerformance(
        prediction_df = prediction_df,
        truth = seu_xen$cell_type_SingleR,
    )

    expect_true(is.list(performance))

    expect_error(evaluatePredictionPerformance(
        prediction_df = prediction_df,
        truth = seu_xen$cell_type_SingleR[1:3],
    ))
})
