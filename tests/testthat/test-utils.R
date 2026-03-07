test_that("Test factorSorted", {
    a <- factorSorted(c("a", "b", "a", "c", "b", "a"), decreasing = FALSE)

    expect_true(is.factor(a))
})
