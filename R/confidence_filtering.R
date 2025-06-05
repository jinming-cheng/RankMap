
#' Optimize Confidence Threshold for Cell Type Prediction
#'
#' Finds the confidence threshold that best balances prediction
#' accuracy and retention rate.
#' Requires ground truth labels for comparison.
#'
#' @param prediction_df A data frame from \code{\link{PredictRankModel}
#'                     (return_confidence = TRUE)}.
#' @param truth A character or factor vector of true labels,
#'              aligned with rows of \code{prediction_df}.
#' @param thresholds Numeric vector of thresholds to test.
#'                   Default is seq(0.1, 0.9, by = 0.05).
#' @param plot Logical. If \code{TRUE}, generates a plot of
#'             accuracy vs. threshold. Default is \code{TRUE}.
#'
#' @return A data frame summarizing accuracy and retained cell count
#'         at each threshold.
#'
#' @examples
#' # Simulated prediction data
#' pred_df <- data.frame(
#'   cell_id = paste0("cell", 1:10),
#'   predicted_cell_type = c("A", "B", "A", "B", "A", "B", "A", "A", "B", "B"),
#'   confidence = c(0.95, 0.87, 0.65, 0.48, 0.92, 0.55, 0.73, 0.33, 0.99, 0.60)
#' )
#'
#' # Ground truth labels
#' truth <- c("A", "B", "A", "B", "B", "B", "A", "A", "B", "B")
#'
#' # Evaluate how accuracy and coverage change with threshold
#' summary_df <- OptimizeConfidenceThreshold(pred_df, truth, plot = TRUE)
#'
#' # View result
#' print(summary_df)
#'
#' @export
OptimizeConfidenceThreshold <- function(prediction_df,
                                        truth,
                                        thresholds = seq(0.1, 0.9, by = 0.05),
                                        plot = TRUE) {

  if (!all(c("predicted_cell_type", "confidence") %in%
           colnames(prediction_df))) {
    stop("prediction_df must contain columns: ",
         "predicted_cell_type and confidence")
  }

  if (length(truth) != nrow(prediction_df)) {
    stop("Length of truth must match number of rows in prediction_df.")
  }

  results <- lapply(thresholds, function(th) {
    pred <- prediction_df$predicted_cell_type
    pred[prediction_df$confidence < th] <- NA
    valid <- !is.na(pred)

    acc <- if (sum(valid) > 0) {
      mean(pred[valid] == truth[valid])
    } else {
      NA
    }

    data.frame(
      threshold = th,
      accuracy = acc,
      retained = sum(valid),
      retained_frac = mean(valid)
    )
  })

  results_df <- do.call(rbind, results)

  if (plot) {
    plot(results_df$threshold, results_df$accuracy, type = "b",
         ylim = c(0, 1), col = "blue", pch = 16,
         xlab = "Confidence Threshold",
         ylab = "Accuracy (Retained Cells Only)")
    graphics::lines(results_df$threshold,
                    results_df$retained_frac,
                    type = "b", col = "darkgray", pch = 1)
    graphics::legend("bottomleft",
                     legend = c("Accuracy", "Fraction Retained"),
                     col = c("blue", "darkgray"),
                     lty = 1, pch = c(16, 1), bty = "n")
  }

  return(results_df)
}



#' Tag Low-Confidence Predictions and Add Status Column
#'
#' Replaces low-confidence cell type predictions with a placeholder
#' label (e.g., "unknown") and adds a confidence status column
#' ("confident" or "uncertain").
#'
#' @param prediction_df A data frame from \code{\link{PredictRankModel}
#'                      (return_confidence = TRUE)}.
#' @param threshold Numeric. Confidence threshold below which predictions
#'                  are flagged. Default is \code{0.5}.
#' @param low_conf_label Character. Replacement for low-confidence predictions.
#'                       Default is \code{"unknown"}.
#' @param keep_confidence Logical. Whether to retain the original
#'                        \code{confidence} column. Default is \code{TRUE}.
#'
#' @return A data frame with columns:
#'         \code{cell_id},
#'         \code{predicted_cell_type},
#'         \code{status},
#'         nd optionally \code{confidence}.
#'
#' @examples
#' # Simulated predictions with confidence
#' pred_df <- data.frame(
#'   cell_id = paste0("cell", 1:5),
#'   predicted_cell_type = c("B", "T", "B", "Myeloid", "T"),
#'   confidence = c(0.92, 0.47, 0.88, 0.33, 0.76)
#' )
#'
#' # Apply threshold of 0.5 to flag low-confidence cells
#' result <- FilterLowConfidenceCells(pred_df, threshold = 0.5)
#'
#' # Show result
#' print(result)
#'
#' # Remove confidence column and use custom label
#' # for low-confidence predictions
#' result2 <- FilterLowConfidenceCells(pred_df, threshold = 0.6,
#'                                     low_conf_label = "low_conf",
#'                                     keep_confidence = FALSE)
#'
#' @export
FilterLowConfidenceCells <- function(prediction_df,
                                     threshold = 0.5,
                                     low_conf_label = "unknown",
                                     keep_confidence = TRUE) {

  required_cols <- c("cell_id", "predicted_cell_type", "confidence")

  msg <- paste0("Input data frame must contain columns: ",
                paste(required_cols, collapse = ", "))
  if (!all(required_cols %in% colnames(prediction_df))) {
    stop(msg)
  }

  df <- prediction_df

  df$status <- ifelse(df$confidence < threshold, "uncertain", "confident")
  df$predicted_cell_type[df$confidence < threshold] <- low_conf_label

  if (!keep_confidence) {
    df <- df[, c("cell_id", "predicted_cell_type", "status")]
  }

  return(df)
}

