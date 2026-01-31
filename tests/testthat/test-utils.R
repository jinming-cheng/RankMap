test_that("Test FactorSorted", {
    a <- FactorSorted(c("a", "b", "a", "c", "b", "a"), decreasing = FALSE)

    expect_true(is.factor(a))
})
