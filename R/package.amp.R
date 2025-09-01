#' @title Load Required R Packages for Data Analysis
#' @description
#' Checks and installs (if missing) various R packages commonly used for microbiome analysis.
#' @note
#' This function will check package availability and attempt to install missing packages from CRAN.
#' @return Invisibly returns a logical vector indicating which packages were successfully loaded
#' @export
#' @examples
#' \dontrun{
#' # Check and load packages
#' package.amp()
#' }
package.amp <- function() {
  # 包列表
  pkg_list <- c(
    "phyloseq", "tidyverse",   "fs",
    "ggthemes", "RColorBrewer", "magrittr", "ggsignif",
    "ggtree", "ggtreeExtra", "ggstar", "MicrobiotaProcess",
    "ggnewscale", "grid", "caret"
  )

  # 检查并安装缺失包
  new_pkgs <- pkg_list[!pkg_list %in% installed.packages()[,"Package"]]
  if(length(new_pkgs)) install.packages(new_pkgs)

  # 加载所有包
  suppressPackageStartupMessages(
    sapply(pkg_list, require, character.only = TRUE)
  )

  invisible(sapply(pkg_list, function(x) x %in% .packages()))
}
