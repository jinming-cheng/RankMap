
#' Sample Cells by Cell Type with Minimum Representation
#'
#' This function samples a specified total number of cells from
#' a metadata table, ensuring that each cell type is represented
#' by at least a minimum number of cells.
#'
#' @param cell_metadata A data frame containing at least two columns:
#'                      \code{cell_id} and \code{cell_type}.
#' @param n_total_cells Total number of cells to sample.
#'                      Default is \code{5000}.
#' @param min_per_type Minimum number of cells to sample from each
#'                     \code{cell_type}. Default is \code{50}.
#'
#' @return A data frame of sampled cells with at least
#'         \code{min_per_type} cells per \code{cell_type}.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @examples
#' meta <- data.frame(
#'   cell_id = paste0("cell", 1:10000),
#'   cell_type = sample(c("T", "B", "Mac"), 10000, replace = TRUE)
#' )
#' set.seed(42)
#' sampled <- SampleCellsByType(meta, n_total_cells = 1000,
#'                              min_per_type = 50)
#'
#' @export
SampleCellsByType <- function(cell_metadata,
                              n_total_cells = 5000,
                              min_per_type = 50) {

  # Validate input
  required_cols <- c("cell_id", "cell_type")
  missing_cols <- setdiff(required_cols, names(cell_metadata))
  msg <- paste0("Input 'cell_metadata' must contain columns: ",
                paste(missing_cols, collapse = ", "))
  if (length(missing_cols) > 0) {
    stop(msg)
  }

  # Ensure each cell type has at least `min_per_type` cells
  initial_cells <- cell_metadata %>%
    dplyr::group_by(.data$cell_type) %>%
    dplyr::slice_sample(n = min_per_type, replace = FALSE) %>%
    dplyr::ungroup()

  n_additional_needed <- n_total_cells - nrow(initial_cells)

  if (n_additional_needed > 0) {
    remaining_data <- dplyr::anti_join(cell_metadata,
                                       initial_cells, by = "cell_id")

    additional_cells <- remaining_data %>%
      dplyr::sample_n(n_additional_needed)

    final_cells <- dplyr::bind_rows(initial_cells, additional_cells)
  } else {
    final_cells <- initial_cells
  }

  return(final_cells)
}

