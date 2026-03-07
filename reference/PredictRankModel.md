# Predict Cell Types, Probabilities, or Confidence Scores

Predicts cell type labels or class probabilities using a trained RankMap
model. Supports both `glmnet` and `cv.glmnet`. Optionally returns a data
frame with predicted cell types and confidence scores.

## Usage

``` r
predictRankModel(
  model,
  new_data,
  lambda = NULL,
  return_probs = FALSE,
  return_confidence = FALSE,
  ...
)
```

## Arguments

- model:

  A fitted model from
  [`trainRankModel`](https://github.com/jinming-cheng/RankMap/reference/TrainRankModel.md).

- new_data:

  A gene expression matrix with genes as rows and cells (or spatial
  spots) as columns. Can be a `Seurat` object, a `SummarizedExperiment`
  object, a dense numeric matrix, or a sparse `dgCMatrix`.

- lambda:

  Optional lambda value. If `NULL`, uses `lambda.min` (if available),
  else `0.01`.

- return_probs:

  Logical. If `TRUE`, returns the full matrix of class probabilities.

- return_confidence:

  Logical. If `TRUE`, returns a data frame with predicted labels and
  confidence scores.

- ...:

  Additional arguments passed to
  [`computeRankedMatrix`](https://github.com/jinming-cheng/RankMap/reference/ComputeRankedMatrix.md).

## Value

A character vector (default), a matrix (if `return_probs = TRUE`), or a
data frame (if `return_confidence = TRUE`).

## Examples

``` r
# Read in single-cell reference data
seu_sc <- readRDS(system.file("extdata", "seu_sc.rds",
    package = "RankMap"
))

# Read in Xenium spatial data
seu_xen <- readRDS(system.file("extdata", "seu_xen.rds",
    package = "RankMap"
))

# Extract normalized expression data
common_genes <- intersect(rownames(seu_sc), rownames(seu_xen))
mat <- extractData(seu_sc)[common_genes, ]
new_mat <- extractData(seu_xen)[common_genes, ]

# Train a model
set.seed(42)
model <- trainRankModel(mat, seu_sc$cell_type)

# Predict cell type
pred <- predictRankModel(model, new_mat)

table(predict = pred, truth = seu_xen$cell_type_SingleR)
#>        truth
#> predict Basal LP Tumor
#>   Basal    48  4     1
#>   LP        2 45     0
#>   Tumor     0  1    49
```
