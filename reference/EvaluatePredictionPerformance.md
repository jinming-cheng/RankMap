# Evaluate Prediction Accuracy and Confusion by Cell Type

Computes overall accuracy, per-class accuracy, and confusion matrix
given predicted vs. true labels.

## Usage

``` r
EvaluatePredictionPerformance(
  prediction_df = NULL,
  truth = NULL,
  low_conf_label = "unknown"
)
```

## Arguments

- prediction_df:

  Output from
  [`PredictRankModel`](https://github.com/jinming-cheng/RankMap/reference/PredictRankModel.md)` (return_confidence = TRUE)`
  or filtered version.

- truth:

  A vector of true labels matching `prediction_df`.

- low_conf_label:

  Character. Label used for low-confidence predictions (e.g.,
  `"unknown"` or `"uncertain"`). These will be excluded from evaluation.
  Set to `NULL` to include all predictions.

## Value

A list with:

- overall_accuracy:

  Proportion of correct predictions

- per_class_accuracy:

  Accuracy per true cell type

- confusion_matrix:

  Contingency table (true × predicted)

## Examples

``` r
# Read in single-cell reference data
seu_sc <- readRDS(system.file("extdata", "seu_sc.rds",
    package = "RankMap"
))
#> Loading required package: SeuratObject
#> Loading required package: sp
#> ‘SeuratObject’ was built under R 4.4.0 but the current version is
#> 4.4.2; it is recomended that you reinstall ‘SeuratObject’ as the ABI
#> for R may have changed
#> 
#> Attaching package: ‘SeuratObject’
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, t

# Read in Xenium spatial data
seu_xen <- readRDS(system.file("extdata", "seu_xen.rds",
    package = "RankMap"
))

# Predict cell type for Xenium data
prediction_df <- RankMap(
    ref_data = seu_sc,
    ref_labels = seu_sc$cell_type,
    new_data = seu_xen
)

performance <- EvaluatePredictionPerformance(
    prediction_df = prediction_df,
    truth = seu_xen$cell_type_SingleR
)
performance
#> $overall_accuracy
#> [1] 0.9466667
#> 
#> $per_class_accuracy
#> Basal    LP Tumor 
#>  0.96  0.90  0.98 
#> 
#> $confusion_matrix
#>        Predicted
#> True    Basal LP Tumor
#>   Basal    48  2     0
#>   LP        4 45     1
#>   Tumor     1  0    49
#> 
```
