#' Compute the arithmetic mean
#'
#' This function computes the arithmetic mean of a numeric variable.
#'
#' @param x a numeric vector
#'
#' @return A `numeric` representing the arithmetic mean of `x`.
#'
#' @export
#'
#' @examples
#' x <- 1:10
#' moyenne(x)

moyenne <- function(x) sum(x) / length(x)
