test_that("Test sampleCellsByType", {
    meta <- data.frame(
        cell_id = paste0("cell", 1:10000),
        cell_type = sample(c("T", "B", "Mac"), 10000, replace = TRUE)
    )

    sampled <- sampleCellsByType(meta,
        n_total_cells = 1000,
        min_per_type = 50
    )

    expect_true(is.data.frame(sampled))

    expect_error(sampleCellsByType(meta[, 1, drop = FALSE],
        n_total_cells = 1000,
        min_per_type = 50
    ))
})
