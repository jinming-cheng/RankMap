# RankMap: Train and Predict Cell Types Using Top-Ranked Genes

A unified function that trains a multinomial classifier on reference
expression data and predicts cell types for a new dataset using top-k
ranked genes. Automatically matches genes between datasets, optionally
performs cross-validation, and returns predictions with confidence
scores or probabilities.

## Usage

``` r
RankMap(
  ref_data = NULL,
  ref_labels = NULL,
  new_data = NULL,
  n_feature_max = 500,
  k = 20,
  alpha = 0.1,
  cv = FALSE,
  nfolds = 5,
  lambda = NULL,
  return_probs = FALSE,
  return_confidence = TRUE,
  threshold = NULL,
  return_model = FALSE,
  ...
)
```

## Arguments

- ref_data:

  Reference gene expression matrix (genes x cells), a `Seurat` object,
  or a `SummarizedExperiment` object.

- ref_labels:

  A character or factor vector of cell type labels for columns of
  `ref_data`.

- new_data:

  New data to annotate. Same format as `ref_data` (matrix, `dgCMatrix`,
  `Seurat` object or `SummarizedExperiment` object).

- n_feature_max:

  Maximum number of genes to use when more than 500 genes are shared.
  Default is `500`.

- k:

  Number of top expressed genes to retain per cell (ranking). Default is
  `20`.

- alpha:

  Elastic net mixing parameter for `glmnet`. Default is `0.1`.

- cv:

  Logical. Whether to use `cv.glmnet` for cross-validation. Default is
  `FALSE`.

- nfolds:

  Number of folds for cross-validation. Default is `5`.

- lambda:

  Optional lambda value for prediction. If `NULL`, uses `lambda.min`
  from CV or defaults to `0.01`.

- return_probs:

  Logical. If `TRUE`, returns full class probability matrix. Default is
  `FALSE`.

- return_confidence:

  Logical. If `TRUE`, returns prediction with confidence score and
  status. Default is `TRUE`.

- threshold:

  Optional numeric threshold. If set and `return_confidence = TRUE`,
  predictions below this confidence are labeled as `"unknown"`.

- return_model:

  Logical. If `TRUE`, returns a list containing both predictions and the
  trained model. Default is `FALSE`.

- ...:

  Additional arguments passed to
  [`ComputeRankedMatrix`](https://github.com/jinming-cheng/RankMap/reference/ComputeRankedMatrix.md).

## Value

A data frame of predictions (by default), or a list with elements
`predictions` and `model` if `return_model = TRUE`.

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

# Predict cell type for spatial data using single-cell data as reference
pred_df <- RankMap(
    ref_data = seu_sc,
    ref_labels = seu_sc$cell_type,
    new_data = seu_xen
)
head(pred_df)
#>   cell_id predicted_cell_type confidence
#> 1    3869               Tumor     0.8829
#> 2    5257               Tumor     0.9612
#> 3    6456               Basal     0.9243
#> 4    8555                  LP     0.8847
#> 5    9243               Basal     0.9911
#> 6   10303               Basal     0.9971
```
