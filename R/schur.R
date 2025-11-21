#' @title Schur Product (Element-wise multiplication)
#'
#' Calculates the element-wise product (Hadamard product) of two matrices.
#'
#' @param a numeric matrix
#' @param b numeric matrix
#' @return numeric matrix resulting from element-wise multiplication
#' @keywords internal
#' @export
schur <- function(a, b) {
  .Call("_EasyMultiOmics_schur", PACKAGE = "EasyMultiOmics", a, b)
}
