#' @title Load Required R Packages for Data Analysis
#'
#' @description
#' The function loads various R libraries commonly used for microbiome analysis, data manipulation,
#' and visualization. It includes packages like `phyloseq`, `ggClusterNet`, `EasyStat`, and more.
#' This function is intended to streamline the initial setup for data analysis by loading multiple
#' libraries in one call, saving time and avoiding repetitive loading of the same libraries in
#' separate scripts.
#' @return
#' This function does not return any value. It is used for loading required libraries.
#' @examples
#' \dontrun{
#' # Load all required packages for microbiome analysis
#' package.amp()
#' }
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njau.edu.cn}
package.amp <- function(){
  #导入R包#-------
  library(phyloseq)
  library(tidyverse)
  library(ggClusterNet)
  library(EasyStat)
  library(fs)
  library(ggthemes)
  library(RColorBrewer)#调色板调用包

  # library(ggplot2)
  # library(dplyr)
  library(magrittr)
  # library(RColorBrewer)
  # devtools::install_github("YuLab-SMU/treeio") # 安装1.7以上版本的才能支持MicrobiotaProcess
  # devtools::install_github("YuLab-SMU/MicrobiotaProcess")
  # BiocManager::install("treeio")
  library(MicrobiotaProcess)
  # library(tibble)
  library(ggsignif)
  library(ggtree)
  library(ggtreeExtra)
  # library(ggplot2)
  library(ggstar)
  library(MicrobiotaProcess)
  library(ggnewscale)
  library(grid)
  library(caret)
}

