#' Convert phyloseq::sample_data to a plain data.frame
#'
#' @param x A phyloseq::sample_data object, or a data.frame.
#' @return A base data.frame with rownames preserved and class stripped to "data.frame".
#' @export
as.data.frame2 <- function(x) {
  # 如果是 sample_data（即使同时也是 data.frame），优先这样处理
  if (inherits(x, "sample_data")) {
    df <- methods::as(x, "data.frame")
  } else if (is.data.frame(x)) {
    df <- x
  } else {
    stop("x must be a phyloseq::sample_data or data.frame")
  }

  # 再保险一次转成 data.frame，并且去掉乱七八糟的 class
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  class(df) <- "data.frame"

  df
}

