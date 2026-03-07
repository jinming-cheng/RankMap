#' Extract Expression Matrix from Seurat, SummarizedExperiment, or Matrix
#'
#' Extracts a gene expression matrix from a \code{Seurat} object
#' (default assay and \code{"data"} slot),
#' a \code{SummarizedExperiment} object (assay named \code{"logcounts"}),
#' or returns the input matrix if it is already a dense or sparse matrix.
#'
#' @param data A \code{Seurat} object, a \code{SummarizedExperiment} object,
#'   a dense numeric matrix, or a sparse matrix of class \code{dgCMatrix}
#'
#' @return A numeric gene expression matrix (dense or sparse),
#'         with genes as rows and cells or spots as columns.
#'
#' @examples
#' seu_sc <- readRDS(system.file("extdata", "seu_sc.rds",
#'     package = "RankMap"
#' ))
#'
#' # From Seurat object:
#' mat <- extractData(seu_sc)
#'
#' @export
extractData <- function(data) {
    if (!inherits(data, c(
        "Seurat", "SummarizedExperiment",
        "matrix", "dgCMatrix"
    ))) {
        msg <- paste0(
            "Input must be a Seurat object, ",
            "a SummarizedExperiment object, ",
            "a numeric matrix, or a dgCMatrix."
        )
        stop(msg)
    }

    if (inherits(data, "Seurat")) {
        return(Seurat::GetAssayData(data))
    }

    if (inherits(data, "SummarizedExperiment")) {
        return(SummarizedExperiment::assay(data, "logcounts"))
    }

    return(data)
}


#' Mask Non-Top-K Genes in Each Cell or Spot
#'
#' This function zeroes out all but the top-k most highly expressed genes
#' in each column of a gene expression matrix.
#' It supports both dense and sparse inputs.
#'
#' @param data A gene expression matrix with genes as rows and cells
#'             (or spatial spots) as columns.
#'             Can be a dense numeric matrix, or a sparse \code{dgCMatrix}.
#' @param k Integer. Number of top expressed genes to retain per column.
#'          Default is \code{20}.
#'
#' @return A dense numeric matrix with the same dimensions as \code{data},
#'         with only the top-k values retained per column.
#' @examples
#' mat <- matrix(runif(1000), nrow = 100)
#' masked <- maskTopKGenes(mat, k = 10)
#'
#' @export
maskTopKGenes <- function(data, k = 20) {
    if (!inherits(data, c("matrix", "dgCMatrix"))) {
        stop("Input 'data' must be a matrix or dgCMatrix.")
    }

    # Use apply() for performance — converts one column at a time to dense
    masked <- apply(data, 2, function(col) {
        top_k_idx <- order(col,
            decreasing = TRUE
        )[seq_len(min(k, length(col)))]
        out <- numeric(length(col))
        out[top_k_idx] <- col[top_k_idx]
        out
    })

    # Preserve orientation
    if (!all(dim(masked) == dim(data))) {
        masked <- t(masked)
    }

    dimnames(masked) <- dimnames(data)
    return(masked)
}


#' Compute Ranked Gene Expression Matrix
#'
#' This function transforms a gene expression matrix by
#' applying top-k filtering, optional rank transformation, binning,
#' scaling, and/or expression weighting.
#'
#' If \code{use_data = TRUE}, the function bypasses all transformation
#' steps and returns the full input expression matrix.
#' This is useful for benchmarking performance
#' of rank-based versus raw log-normalized expression models.
#'
#' @param data A gene expression matrix with genes as rows and cells
#'             (or spatial spots) as columns.
#'             Can be a dense numeric matrix, or a sparse \code{dgCMatrix}.
#' @param weight_by_expr Logical. Whether to weight the ranks by
#'                       log-transformed expression values.
#'                       Default is \code{TRUE}.
#' @param rank_zeros Logical. Whether to include zero values in the ranking.
#'                   Default is \code{FALSE}.
#' @param bin_rank Logical. Whether to discretize ranks into bins.
#'                          Default is \code{TRUE}.
#' @param scale_rank Logical. Whether to z-score normalize the ranks
#'                   across columns. Default is \code{TRUE}.
#' @param k Integer. Number of top expressed genes to retain per column.
#'                   Default is \code{20}.
#' @param use_data Logical. If \code{TRUE}, returns the full input
#'                 expression matrix (unfiltered and unranked).
#'                 Default is \code{FALSE}.
#'
#' @return A numeric matrix of ranked expression values
#'         (or the original expression matrix if \code{use_data = TRUE}).
#'         Output is always returned as a dense \code{matrix}.
#'
#' @seealso \code{\link{maskTopKGenes}},
#'          \code{\link{trainRankModel}},
#'          \code{\link{predictRankModel}}
#'
#' @examples
#' mat <- matrix(runif(1000), nrow = 100)
#' ranked <- computeRankedMatrix(mat, k = 10)
#' raw_expr <- computeRankedMatrix(mat, use_data = TRUE)
#'
#' @export
computeRankedMatrix <- function(
    data,
    weight_by_expr = TRUE,
    rank_zeros = FALSE,
    bin_rank = TRUE,
    scale_rank = TRUE,
    k = 20,
    use_data = FALSE) {
    if (use_data) {
        return(data)
    }

    data <- maskTopKGenes(data, k = k)

    # Replace zeros with NA if rank_zeros is FALSE
    if (!rank_zeros) {
        data[data == 0] <- NA
    }

    # Rank within each column (cell/spot)
    ranked <- apply(data, 2, rank, ties.method = "average", na.last = "keep")

    # Optional binning of ranks
    if (bin_rank) {
        ranked <- apply(ranked, 1, function(row) {
            cut(row, breaks = seq(0, k, length.out = 5), labels = FALSE)
        })
        ranked <- t(ranked)
    }

    # Replace NAs back to zero
    if (!rank_zeros) {
        ranked[is.na(ranked)] <- 0
    }

    # Optional weighting by log-transformed expression
    if (weight_by_expr) {
        if (!rank_zeros) data[is.na(data)] <- 0
        ranked <- ranked * log1p(round(data, 2))
    }

    # Optional z-score normalization across columns
    if (scale_rank) {
        ranked <- scale(ranked)
    }

    return(ranked)
}
