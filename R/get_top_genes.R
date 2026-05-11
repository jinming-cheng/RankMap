#' Get top expressed genes for each cell type
#'
#' Identify the top \code{n} highly expressed genes for each cell type
#' based on average expression across cells. This function is useful for
#' summarizing representative genes from a reference dataset and can be
#' optionally used to restrict features for downstream analysis.
#'
#' @param data A gene expression object. Supported inputs include
#' \code{Seurat}, \code{SingleCellExperiment}, \code{SpatialExperiment},
#' or a matrix-like object. Expression values are extracted using
#' \code{extractData()}.
#'
#' @param cell_type A vector of cell type labels corresponding to columns
#' of \code{data}.
#'
#' @param n Integer. Number of top expressed genes to return for each
#' cell type. Default is \code{50}.
#'
#' @param genes Optional character vector of gene names to restrict the
#' analysis. Only genes present in \code{data} will be used.
#'
#' @param return_list Logical. If \code{TRUE}, return a list of top genes
#' for each cell type. If \code{FALSE}, return a unique vector of genes
#' across all cell types. Default is \code{FALSE}.
#'
#' @return
#' If \code{return_list = TRUE}, a named list where each element contains
#' the top \code{n} genes for a cell type. Otherwise, a character vector
#' of unique top genes across all cell types.
#'
#' @details
#' For each cell type, genes are ranked by their average expression
#' (computed using \code{Matrix::rowMeans}) across cells of that type.
#' Cells with missing (\code{NA}) cell type labels are removed prior to
#' computation.
#'
#' @examples
#' # Read in single-cell reference data
#' seu_sc <- readRDS(system.file("extdata", "seu_sc.rds",
#'     package = "RankMap"
#' ))
#'
#' # Get top genes across all cell types
#' top_genes <- topCellTypeGenes(
#'     data = seu_sc,
#'     cell_type = seu_sc$cell_type,
#'     n = 50
#' )
#'
#' # Get top genes per cell type
#' top_list <- topCellTypeGenes(
#'     data = seu_sc,
#'     cell_type = seu_sc$cell_type,
#'     n = 50,
#'     return_list = TRUE
#' )
#'
#' @export
topCellTypeGenes <- function(
    data,
    cell_type,
    n = 50,
    genes = NULL,
    return_list = FALSE
) {
    mat <- extractData(data)

    if (length(cell_type) != ncol(mat)) {
        stop(
            "Length of 'cell_type' must match ",
            "the number of columns in 'data'."
        )
    }

    if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 1) {
        stop("'n' must be a positive integer.")
    }
    n <- as.integer(n)

    if (!is.null(genes)) {
        genes <- intersect(genes, rownames(mat))
        if (length(genes) == 0) {
            stop("None of the input 'genes' are found in the data.")
        }
        mat <- mat[genes, , drop = FALSE]
    }

    keep <- !is.na(cell_type)
    if (!all(keep)) {
        mat <- mat[, keep, drop = FALSE]
        cell_type <- cell_type[keep]
        warning("Cells with NA cell_type were removed.")
    }

    ct <- factor(cell_type)
    idx <- split(seq_len(ncol(mat)), ct)

    top_list <- lapply(idx, function(j) {
        avg_expr <- Matrix::rowMeans(mat[, j, drop = FALSE])
        avg_expr <- sort(avg_expr, decreasing = TRUE)
        names(avg_expr)[seq_len(min(n, length(avg_expr)))]
    })

    if (return_list) {
        return(top_list)
    }

    unique(unlist(top_list, use.names = FALSE))
}
