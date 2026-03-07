# Tag Low-Confidence Predictions and Add Status Column

Replaces low-confidence cell type predictions with a placeholder label
(e.g., "unknown") and adds a confidence status column ("confident" or
"uncertain").

## Usage

``` r
filterLowConfidenceCells(
  prediction_df,
  threshold = 0.5,
  low_conf_label = "unknown",
  keep_confidence = TRUE
)
```

## Arguments

- prediction_df:

  A data frame from
  [`predictRankModel`](https://github.com/jinming-cheng/RankMap/reference/PredictRankModel.md)` (return_confidence = TRUE)`.

- threshold:

  Numeric. Confidence threshold below which predictions are flagged.
  Default is `0.5`.

- low_conf_label:

  Character. Replacement for low-confidence predictions. Default is
  `"unknown"`.

- keep_confidence:

  Logical. Whether to retain the original `confidence` column. Default
  is `TRUE`.

## Value

A data frame with columns: `cell_id`, `predicted_cell_type`, `status`,
nd optionally `confidence`.

## Examples

``` r
# Simulated predictions with confidence
pred_df <- data.frame(
    cell_id = paste0("cell", 1:5),
    predicted_cell_type = c("B", "T", "B", "Myeloid", "T"),
    confidence = c(0.92, 0.47, 0.88, 0.33, 0.76)
)

# Apply threshold of 0.5 to flag low-confidence cells
result <- filterLowConfidenceCells(pred_df, threshold = 0.5)

# Show result
print(result)
#>   cell_id predicted_cell_type confidence    status
#> 1   cell1                   B       0.92 confident
#> 2   cell2             unknown       0.47 uncertain
#> 3   cell3                   B       0.88 confident
#> 4   cell4             unknown       0.33 uncertain
#> 5   cell5                   T       0.76 confident

# Remove confidence column and use custom label
# for low-confidence predictions
result2 <- filterLowConfidenceCells(pred_df,
    threshold = 0.6,
    low_conf_label = "low_conf",
    keep_confidence = FALSE
)
```
