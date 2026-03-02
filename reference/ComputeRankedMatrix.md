# Compute Ranked Gene Expression Matrix

This function transforms a gene expression matrix by applying top-k
filtering, optional rank transformation, binning, scaling, and/or
expression weighting.

## Usage

``` r
ComputeRankedMatrix(
  data,
  weight_by_expr = TRUE,
  rank_zeros = FALSE,
  bin_rank = TRUE,
  scale_rank = TRUE,
  k = 20,
  use_data = FALSE
)
```

## Arguments

- data:

  A gene expression matrix with genes as rows and cells (or spatial
  spots) as columns. Can be a dense numeric matrix, or a sparse
  `dgCMatrix`.

- weight_by_expr:

  Logical. Whether to weight the ranks by log-transformed expression
  values. Default is `TRUE`.

- rank_zeros:

  Logical. Whether to include zero values in the ranking. Default is
  `FALSE`.

- bin_rank:

  Logical. Whether to discretize ranks into bins. Default is `TRUE`.

- scale_rank:

  Logical. Whether to z-score normalize the ranks across columns.
  Default is `TRUE`.

- k:

  Integer. Number of top expressed genes to retain per column. Default
  is `20`.

- use_data:

  Logical. If `TRUE`, returns the full input expression matrix
  (unfiltered and unranked). Default is `FALSE`.

## Value

A numeric matrix of ranked expression values (or the original expression
matrix if `use_data = TRUE`). Output is always returned as a dense
`matrix`.

## Details

If `use_data = TRUE`, the function bypasses all transformation steps and
returns the full input expression matrix. This is useful for
benchmarking performance of rank-based versus raw log-normalized
expression models.

## See also

[`MaskTopKGenes`](https://github.com/jinming-cheng/RankMap/reference/MaskTopKGenes.md),
[`TrainRankModel`](https://github.com/jinming-cheng/RankMap/reference/TrainRankModel.md),
[`PredictRankModel`](https://github.com/jinming-cheng/RankMap/reference/PredictRankModel.md)

## Examples

``` r
mat <- matrix(runif(1000), nrow = 100)
ranked <- ComputeRankedMatrix(mat, k = 10)
raw_expr <- ComputeRankedMatrix(mat, use_data = TRUE)
```
