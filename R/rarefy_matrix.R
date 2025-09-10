#' Rarefy a count matrix
#'
#' Rarefy the OTU/feature matrix to depth N.
#'
#' @param mat numeric matrix of counts
#' @param N integer, rarefaction depth
#' @return numeric matrix of rarefied counts
#' @export
rarefy_matrix <- function(mat, N) {
  .Call("_EasyMultiOmics_rarefy_matrix", PACKAGE = "EasyMultiOmics", mat, N)
}
