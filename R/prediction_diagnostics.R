#' Evaluate Prediction Accuracy and Confusion by Cell Type
#'
#' Computes overall accuracy, per-class accuracy, and confusion matrix
#' given predicted vs. true labels.
#'
#' @param prediction_df Output from \code{\link{PredictRankModel}
#'                      (return_confidence = TRUE)} or filtered version.
#' @param truth A vector of true labels matching \code{prediction_df}.
#'
#' @return A list with:
#'   \item{overall_accuracy}{Proportion of correct predictions}
#'   \item{per_class_accuracy}{Accuracy per true cell type}
#'   \item{confusion_matrix}{Contingency table (true × predicted)}
#'
#' @examples
#' # Read in single-cell reference data
#' seu_sc <- readRDS(system.file("extdata", "seu_sc.rds",
#'     package = "RankMap"
#' ))
#'
#' # Read in Xenium spatial data
#' seu_xen <- readRDS(system.file("extdata", "seu_xen.rds",
#'     package = "RankMap"
#' ))
#'
#' # Predict cell type for Xenium data
#' prediction_df <- RankMap(
#'     ref_data = seu_sc,
#'     ref_labels = seu_sc$cell_type,
#'     new_data = seu_xen
#' )
#'
#' performance <- EvaluatePredictionPerformance(
#'     prediction_df = prediction_df,
#'     truth = seu_xen$cell_type_SingleR
#' )
#' performance
#'
#' @export
EvaluatePredictionPerformance <- function(prediction_df = NULL, truth = NULL) {
    if (length(truth) != nrow(prediction_df)) {
        stop("Truth labels must match the number of predicted rows.")
    }

    pred <- prediction_df$predicted_cell_type
    valid <- !is.na(pred) & pred != "unknown" & pred != "uncertain"

    overall_accuracy <- mean(pred[valid] == truth[valid])

    per_class_accuracy <- tapply(
        pred[valid] == truth[valid],
        truth[valid],
        mean
    )

    confusion <- table(True = truth[valid], Predicted = pred[valid])

    list(
        overall_accuracy = overall_accuracy,
        per_class_accuracy = per_class_accuracy,
        confusion_matrix = confusion
    )
}
