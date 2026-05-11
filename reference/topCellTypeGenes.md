# Get top expressed genes for each cell type

Identify the top `n` highly expressed genes for each cell type based on
average expression across cells. This function is useful for summarizing
representative genes from a reference dataset and can be optionally used
to restrict features for downstream analysis.

## Usage

``` r
topCellTypeGenes(data, cell_type, n = 50, genes = NULL, return_list = FALSE)
```

## Arguments

- data:

  A gene expression object. Supported inputs include `Seurat`,
  `SingleCellExperiment`, `SpatialExperiment`, or a matrix-like object.
  Expression values are extracted using
  [`extractData()`](https://github.com/jinming-cheng/RankMap/reference/ExtractData.md).

- cell_type:

  A vector of cell type labels corresponding to columns of `data`.

- n:

  Integer. Number of top expressed genes to return for each cell type.
  Default is `50`.

- genes:

  Optional character vector of gene names to restrict the analysis. Only
  genes present in `data` will be used.

- return_list:

  Logical. If `TRUE`, return a list of top genes for each cell type. If
  `FALSE`, return a unique vector of genes across all cell types.
  Default is `FALSE`.

## Value

If `return_list = TRUE`, a named list where each element contains the
top `n` genes for a cell type. Otherwise, a character vector of unique
top genes across all cell types.

## Details

For each cell type, genes are ranked by their average expression
(computed using
[`Matrix::rowMeans`](https://rdrr.io/pkg/Matrix/man/colSums-methods.html))
across cells of that type. Cells with missing (`NA`) cell type labels
are removed prior to computation.

## Examples

``` r
# Read in single-cell reference data
seu_sc <- readRDS(system.file("extdata", "seu_sc.rds",
    package = "RankMap"
))

# Get top genes across all cell types
top_genes <- topCellTypeGenes(
    data = seu_sc,
    cell_type = seu_sc$cell_type,
    n = 50
)

# Get top genes per cell type
top_list <- topCellTypeGenes(
    data = seu_sc,
    cell_type = seu_sc$cell_type,
    n = 50,
    return_list = TRUE
)
```
