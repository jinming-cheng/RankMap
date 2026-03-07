# Train Multinomial Rank-Based Model

This function trains a multinomial logistic regression model using
ranked gene expression. Optionally performs cross-validation to select
the optimal regularization parameter.

## Usage

``` r
trainRankModel(
  data = NULL,
  labels = NULL,
  alpha = 0.1,
  cv = FALSE,
  nfolds = 5,
  ...
)
```

## Arguments

- data:

  A gene expression matrix with genes as rows and cells (or spatial
  spots) as columns. Can be a `Seurat` object, a `SummarizedExperiment`
  object, a dense numeric matrix, or a sparse `dgCMatrix`.

- labels:

  A vector of cell type labels (one per column of `data`).

- alpha:

  Elastic net mixing parameter. 1 = LASSO, 0 = Ridge. Default is `0.1`.

- cv:

  Logical. If `TRUE`, performs cross-validation using `cv.glmnet`.
  Default is `FALSE`.

- nfolds:

  Number of folds for cross-validation (only used if `cv = TRUE`).
  Default is `5`.

- ...:

  Additional arguments passed to
  [`computeRankedMatrix`](https://github.com/jinming-cheng/RankMap/reference/ComputeRankedMatrix.md).

## Value

A fitted `glmnet` or `cv.glmnet` model object.

## See also

[`computeRankedMatrix`](https://github.com/jinming-cheng/RankMap/reference/ComputeRankedMatrix.md),
[`predictRankModel`](https://github.com/jinming-cheng/RankMap/reference/PredictRankModel.md)

## Examples

``` r
# Read in single-cell reference data
seu_sc <- readRDS(system.file("extdata", "seu_sc.rds",
    package = "RankMap"
))

# Extract normalized expression data
mat <- extractData(seu_sc)

# Train a model
set.seed(42)
model <- trainRankModel(mat, seu_sc$cell_type)

# Train a model with cross-validation
model <- trainRankModel(
    data = mat,
    labels = seu_sc$cell_type,
    cv = TRUE,
    nfolds = 3
)
```
