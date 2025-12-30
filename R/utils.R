
#' Order factor levels by frequency
#'
#' Takes a vector, converts it to a factor, and reorders its levels
#' by their frequency of occurrence (most to least frequent by default).
#'
#' @param x A vector (character or factor) to be reordered by level frequency.
#' @param decreasing Logical. If TRUE (default), levels are ordered from most
#'   to least frequent; if FALSE, from least to most frequent.
#'
#' @return A factor with levels ordered by frequency.
#'
#' @details
#' Frequencies are computed with \code{table(x)}, which ignores \code{NA} values
#' by default. \code{NA} entries in \code{x} are preserved in the output but are
#' not considered a level.
#'
#' @examples
#' FactorSorted(c("a","b","a","c","b","a"))
#' FactorSorted(c("a","b","a","c","b","a"), decreasing = FALSE)
#' FactorSorted(factor(c("x","y","x", NA)))
#'
#' @export
FactorSorted <- function(x, decreasing = TRUE){
  x <- factor(x)
  x <- droplevels(x)
  ordered_levels <- sort(table(x), decreasing = decreasing)
  x <- factor(x, levels = names(ordered_levels))
  x
}
