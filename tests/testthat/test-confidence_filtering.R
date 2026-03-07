test_that("Test optimizeConfidenceThreshold and filterLowConfidenceCells", {
    seu_sc <- readRDS(system.file("extdata", "seu_sc.rds",
        package = "RankMap"
    ))

    seu_xen <- readRDS(system.file("extdata", "seu_xen.rds",
        package = "RankMap"
    ))

    pred_df <- RankMap(
        ref_data = seu_sc,
        ref_labels = seu_sc$cell_type,
        new_data = seu_xen
    )

    truth <- seu_xen$cell_type_SingleR
    summary_df <- optimizeConfidenceThreshold(pred_df, truth, plot = FALSE)

    expect_true(is.data.frame(summary_df))

    expect_silent(optimizeConfidenceThreshold(pred_df, truth, plot = TRUE))

    expect_error(optimizeConfidenceThreshold(pred_df, truth[1:3], plot = FALSE))

    expect_error(optimizeConfidenceThreshold(
        pred_df[, 1, drop = FALSE],
        truth[1:3],
        plot = FALSE
    ))


    result <- filterLowConfidenceCells(pred_df)

    expect_true(is.data.frame(result))

    expect_error(filterLowConfidenceCells(pred_df[, 1, drop = FALSE]))
})
