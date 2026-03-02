# Mask Non-Top-K Genes in Each Cell or Spot

This function zeroes out all but the top-k most highly expressed genes
in each column of a gene expression matrix. It supports both dense and
sparse inputs.

## Usage

``` r
MaskTopKGenes(data, k = 20)
```

## Arguments

- data:

  A gene expression matrix with genes as rows and cells (or spatial
  spots) as columns. Can be a dense numeric matrix, or a sparse
  `dgCMatrix`.

- k:

  Integer. Number of top expressed genes to retain per column. Default
  is `20`.

## Value

A dense numeric matrix with the same dimensions as `data`, with only the
top-k values retained per column.

## Examples

``` r
mat <- matrix(runif(1000), nrow = 100)
masked <- MaskTopKGenes(mat, k = 10)
```
