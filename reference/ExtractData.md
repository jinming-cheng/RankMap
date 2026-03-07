# Extract Expression Matrix from Seurat, SummarizedExperiment, or Matrix

Extracts a gene expression matrix from a `Seurat` object (default assay
and `"data"` slot), a `SummarizedExperiment` object (assay named
`"logcounts"`), or returns the input matrix if it is already a dense or
sparse matrix.

## Usage

``` r
extractData(data)
```

## Arguments

- data:

  A `Seurat` object, a `SummarizedExperiment` object, a dense numeric
  matrix, or a sparse matrix of class `dgCMatrix`

## Value

A numeric gene expression matrix (dense or sparse), with genes as rows
and cells or spots as columns.

## Examples

``` r
seu_sc <- readRDS(system.file("extdata", "seu_sc.rds",
    package = "RankMap"
))

# From Seurat object:
mat <- extractData(seu_sc)
```
