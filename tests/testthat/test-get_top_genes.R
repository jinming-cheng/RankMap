test_that("Test topCellTypeGenes", {
    seu_sc <- readRDS(system.file("extdata", "seu_sc.rds",
        package = "RankMap"
    ))

    top_list <- topCellTypeGenes(
        data = seu_sc,
        cell_type = seu_sc$cell_type,
        n = 50,
        return_list = TRUE
    )
    expect_true(is.list(top_list))
    expect_equal(length(top_list[[1]]), 50)

    top_genes <- topCellTypeGenes(
        data = seu_sc,
        cell_type = seu_sc$cell_type,
        genes = rownames(seu_sc),
        n = 50
    )
    expect_true(is.vector(top_genes))

    expect_error(topCellTypeGenes(
        data = seu_sc,
        cell_type = seu_sc$cell_type[1],
        n = 50,
        return_list = TRUE
    ))

    expect_error(topCellTypeGenes(
        data = seu_sc,
        cell_type = seu_sc$cell_type,
        n = "not_a_number",
        return_list = TRUE
    ))

    expect_error(topCellTypeGenes(
        data = seu_sc,
        cell_type = seu_sc$cell_type,
        n = 50,
        genes = c("non_gene_1", "non_gene_2")
    ))

    expect_warning(topCellTypeGenes(
        data = seu_sc,
        cell_type = c(seu_sc$cell_type[-1], NA),
        n = 50
    ))
})
