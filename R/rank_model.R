#' Train Multinomial Rank-Based Model
#'
#' This function trains a multinomial logistic regression model
#' using ranked gene expression. Optionally performs cross-validation
#' to select the optimal regularization parameter.
#'
#' @importFrom glmnet glmnet cv.glmnet
#' @param data A gene expression matrix with genes as rows and cells
#'             (or spatial spots) as columns.
#'             Can be a \code{Seurat} object,
#'             a \code{SummarizedExperiment} object,
#'             a dense numeric matrix, or a sparse \code{dgCMatrix}.
#' @param labels A vector of cell type labels (one per column of \code{data}).
#' @param alpha Elastic net mixing parameter. 1 = LASSO, 0 = Ridge.
#'              Default is \code{0.1}.
#' @param cv Logical. If \code{TRUE}, performs cross-validation using
#'           \code{cv.glmnet}. Default is \code{FALSE}.
#' @param nfolds Number of folds for cross-validation
#'               (only used if \code{cv = TRUE}). Default is \code{5}.
#' @param ... Additional arguments passed to \code{\link{ComputeRankedMatrix}}.
#'
#' @return A fitted \code{glmnet} or \code{cv.glmnet} model object.
#'
#' @seealso \code{\link{ComputeRankedMatrix}}, \code{\link{PredictRankModel}}
#'
#' @examples
#' # Read in single-cell reference data
#' seu_sc <- readRDS(system.file("extdata", "seu_sc.rds",
#'     package = "RankMap"
#' ))
#'
#' # Extract normalized expression data
#' mat <- ExtractData(seu_sc)
#'
#' # Train a model
#' set.seed(42)
#' model <- TrainRankModel(mat, seu_sc$cell_type)
#'
#' # Train a model with cross-validation
#' model <- TrainRankModel(
#'     data = mat,
#'     labels = seu_sc$cell_type,
#'     cv = TRUE,
#'     nfolds = 3
#' )
#'
#' @export
TrainRankModel <- function(
    data = NULL,
    labels = NULL,
    alpha = 0.1,
    cv = FALSE,
    nfolds = 5,
    ...) {
    data <- ExtractData(data)

    # Compute ranked expression
    rank_expr <- ComputeRankedMatrix(data, ...)

    # Transpose to cells (rows) × genes (columns)
    rank_expr <- t(rank_expr)

    labels <- as.factor(labels)

    if (cv) {
        model <- glmnet::cv.glmnet(
            x = rank_expr,
            y = labels,
            family = "multinomial",
            alpha = alpha,
            nfolds = nfolds
        )
    } else {
        model <- glmnet::glmnet(
            x = rank_expr,
            y = labels,
            family = "multinomial",
            alpha = alpha
        )
    }

    return(model)
}


#' Predict Cell Types, Probabilities, or Confidence Scores
#'
#' Predicts cell type labels or class probabilities using a trained
#' RankMap model. Supports both \code{glmnet} and \code{cv.glmnet}.
#' Optionally returns a data frame with predicted cell types
#' and confidence scores.
#'
#' @importFrom stats predict
#' @param model A fitted model from \code{\link{TrainRankModel}}.
#' @param new_data A gene expression matrix with genes as rows and cells
#'                 (or spatial spots) as columns.
#'                 Can be a \code{Seurat} object,
#'                 a \code{SummarizedExperiment} object,
#'                 a dense numeric matrix, or a sparse \code{dgCMatrix}.
#' @param lambda Optional lambda value. If \code{NULL},
#'               uses \code{lambda.min} (if available), else \code{0.01}.
#' @param return_probs Logical. If \code{TRUE}, returns the
#'                     full matrix of class probabilities.
#' @param return_confidence Logical. If \code{TRUE}, returns a data frame
#'                          with predicted labels and confidence scores.
#' @param ... Additional arguments passed to
#'            \code{\link{ComputeRankedMatrix}}.
#'
#' @return A character vector (default),
#'         a matrix (if \code{return_probs = TRUE}),
#'         or a data frame (if \code{return_confidence = TRUE}).
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
#' # Extract normalized expression data
#' common_genes <- intersect(rownames(seu_sc), rownames(seu_xen))
#' mat <- ExtractData(seu_sc)[common_genes, ]
#' new_mat <- ExtractData(seu_xen)[common_genes, ]
#'
#' # Train a model
#' set.seed(42)
#' model <- TrainRankModel(mat, seu_sc$cell_type)
#'
#' # Predict cell type
#' pred <- PredictRankModel(model, new_mat)
#'
#' table(predict = pred, truth = seu_xen$cell_type_SingleR)
#'
#' @export
PredictRankModel <- function(
    model,
    new_data,
    lambda = NULL,
    return_probs = FALSE,
    return_confidence = FALSE,
    ...) {
    new_data <- ExtractData(new_data)

    msg <- paste0(
        "Choose either return_probs = TRUE or ",
        "return_confidence = TRUE, not both."
    )
    if (return_probs && return_confidence) {
        stop(msg)
    }

    rank_expr <- ComputeRankedMatrix(new_data, ...)
    rank_expr <- t(rank_expr)

    lambda_to_use <- lambda
    if (is.null(lambda)) {
        if (inherits(model, "cv.glmnet") && !is.null(model$lambda.min)) {
            lambda_to_use <- model$lambda.min
        } else {
            lambda_to_use <- 0.01
        }
    }

    pred_type <- if (return_probs || return_confidence) "response" else "class"
    pred <- predict(model,
        newx = rank_expr,
        s = lambda_to_use, type = pred_type
    )

    # Handle label prediction
    if (!return_probs && !return_confidence) {
        return(as.vector(pred))
    }

    # Handle probability or confidence output
    if (length(dim(pred)) == 3) {
        prob_mat <- pred[, , 1]
    } else {
        prob_mat <- pred
    }

    if (return_probs) {
        return(prob_mat)
    }

    # Return confidence + label
    max_idx <- max.col(prob_mat, ties.method = "first")
    labels <- colnames(prob_mat)[max_idx]
    confidence <- prob_mat[cbind(seq_len(nrow(prob_mat)), max_idx)]

    out <- data.frame(
        cell_id = colnames(new_data),
        predicted_cell_type = labels,
        confidence = round(confidence, 4),
        stringsAsFactors = FALSE
    )

    return(out)
}



#' RankMap: Train and Predict Cell Types Using Top-Ranked Genes
#'
#' A unified function that trains a multinomial classifier on
#' reference expression data and predicts cell types for a new dataset
#' using top-k ranked genes. Automatically matches genes between datasets,
#' optionally performs cross-validation, and returns predictions
#' with confidence scores or probabilities.
#'
#' @param ref_data Reference gene expression matrix (genes x cells),
#'                 a \code{Seurat} object,
#'                 or a \code{SummarizedExperiment} object.
#' @param ref_labels A character or factor vector of cell type labels for
#'               columns of \code{ref_data}.
#' @param new_data New data to annotate. Same format as \code{ref_data}
#'                 (matrix, \code{dgCMatrix}, \code{Seurat} object
#'                 or \code{SummarizedExperiment} object).
#' @param n_feature_max Maximum number of genes to use when more than
#'                      500 genes are shared. Default is \code{500}.
#' @param k Number of top expressed genes to retain per cell (ranking).
#'          Default is \code{20}.
#' @param alpha Elastic net mixing parameter for \code{glmnet}.
#'              Default is \code{0.1}.
#' @param cv Logical. Whether to use \code{cv.glmnet} for
#'           cross-validation. Default is \code{FALSE}.
#' @param nfolds Number of folds for cross-validation. Default is \code{5}.
#' @param lambda Optional lambda value for prediction. If \code{NULL},
#'               uses \code{lambda.min} from CV or defaults to \code{0.01}.
#' @param return_probs Logical. If \code{TRUE},
#'                     returns full class probability matrix.
#'                     Default is \code{FALSE}.
#' @param return_confidence Logical. If \code{TRUE}, returns prediction
#'                          with confidence score and status.
#'                          Default is \code{TRUE}.
#' @param threshold Optional numeric threshold.
#'                  If set and \code{return_confidence = TRUE},
#'                  predictions below this confidence are labeled
#'                  as \code{"unknown"}.
#' @param return_model Logical. If \code{TRUE}, returns a list containing
#'                     both predictions and the trained model.
#'                     Default is \code{FALSE}.
#' @param ... Additional arguments passed to
#'            \code{\link{ComputeRankedMatrix}}.
#'
#' @return A data frame of predictions (by default),
#'         or a list with elements \code{predictions} and
#'         \code{model} if \code{return_model = TRUE}.
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
#' # Predict cell type for spatial data using single-cell data as reference
#' pred_df <- RankMap(
#'     ref_data = seu_sc,
#'     ref_labels = seu_sc$cell_type,
#'     new_data = seu_xen
#' )
#' head(pred_df)
#'
#' @export
RankMap <- function(
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
    ...) {
    # Extract matrices
    ref_mat <- ExtractData(ref_data)
    test_mat <- ExtractData(new_data)

    # Find common genes
    common_genes <- intersect(rownames(ref_mat), rownames(test_mat))
    if (length(common_genes) < 50) {
        stop("Too few common genes between reference and new data (< 50).")
    }

    # Reduce to variable genes if too many
    if (length(common_genes) > n_feature_max) {
        message(
            "More than ", n_feature_max,
            " common genes found. Using top ", n_feature_max,
            " most variable genes for training and prediction."
        )
        ref_subset <- ref_mat[common_genes, , drop = FALSE]
        gene_vars <- matrixStats::rowVars(as.matrix(ref_subset))
        top_idx <- order(gene_vars, decreasing = TRUE)[seq_len(n_feature_max)]
        common_genes <- rownames(ref_subset)[top_idx]
    }


    # Subset both matrices to common genes
    ref_mat <- ref_mat[common_genes, , drop = FALSE]
    test_mat <- test_mat[common_genes, , drop = FALSE]

    # Train model
    model <- TrainRankModel(
        data = ref_mat,
        labels = ref_labels,
        k = k,
        alpha = alpha,
        cv = cv,
        nfolds = nfolds,
        ...
    )

    # Predict
    pred_df <- PredictRankModel(
        model = model,
        new_data = test_mat,
        return_confidence = return_confidence,
        return_probs = return_probs,
        k = k,
        lambda = lambda,
        ...
    )

    # Apply confidence threshold
    if (return_confidence && !is.null(threshold)) {
        pred_df <- FilterLowConfidenceCells(pred_df, threshold = threshold)
    }

    if (return_model) {
        return(list(predictions = pred_df, model = model))
    } else {
        return(pred_df)
    }
}
