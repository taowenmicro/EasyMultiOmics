#' @title Calculate and Visualize Bray-Curtis Distances Between Groups
#'
#' @description
#' The `distance.micro` function computes Bray-Curtis distances between groups of microbial community samples.
#' It generates pairwise group comparisons and visualizes the results using boxplots and bar plots with statistical annotations.
#'
#' @param ps A `phyloseq` object containing microbial community data.
#' @param group A character string specifying the grouping variable in the metadata (default: "Group").
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{p1_1}: A ggplot object showing boxplots of Bray-Curtis distances.
#'   \item \code{p1_2}: A ggplot object showing bar plots of Bray-Curtis distances.
#'   \item \code{distance.data}: A data frame with pairwise distances and corresponding group comparisons.
#' }
#'
#' @details
#' This function calculates the Bray-Curtis distances between all samples in the dataset using the `vegan` package.
#' For each pair of groups specified by the `group` variable, it extracts the distances between samples
#' and organizes the results into a data frame. The function then visualizes the distances using boxplots
#' and bar plots with statistical annotations added through the `EasyStat` package.
#'
#' @examples
#' # Example usage:
#' library(phyloseq)
#' library(EasyStat)
#' library(vegan)
#'
#' # Assuming `ps` is a phyloseq object with grouping variable "Treatment"
#' result <- distance.micro(ps = ps, group = "Treatment")
#'
#' # Access and display the plots
#' result[[1]]  # Boxplot
#' result[[2]]  # Bar plot
#' result[[3]]  # Combined plot
#'
#' # Access the distance data
#' head(result$distance.data)
#'
#' @export

distance.micro = function(
    ps =NULL,
    group = "Group"
  ){

  map = as.data.frame(sample_data(ps))
  gro = map[,group] %>% unique()
  colnames(gro) = "group"
  conbgroup = combn(gro$group,2)
  # 计算包括终点均值的所有样品bray距离
  bray_curtis = vegan::vegdist(vegan_otu(ps), method = "bray")
  bray_curtis = as.matrix(bray_curtis)

  for (i in 1:dim(conbgroup)[2]) {
    a = conbgroup[,i]
    map = as.data.frame(sample_data(ps))
    head(map)

    chose1 = map[as.matrix(map[,group]) %>% as.vector() == a[1],] %>% row.names()
    chose2 = map[as.matrix(map[,group]) %>% as.vector() == a[2],] %>% row.names()

    dat = data.frame(group = paste(a[1],a[2],sep = "_VS_"), Distance =bray_curtis[chose1,chose2] %>% as.dist() %>% as.vector() )
    head(dat)

    if (i == 1) {
      table = dat
    }

    if (i != 1) {
      table = rbind(table,dat)
    }
  }

  head(table)
  table$id = 1:dim(table)[1]
  data <- table %>% dplyr::select(id,everything())


  result = EasyStat::MuiKwWlx(data = data,num = c(3))
  result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = c(3),result = result,sig_show ="abc",ncol = 1 )
  p1_1 = result1[[1]]
  p1_1
  res = EasyStat::FacetMuiPlotresultBar(data = data,num = c(3),result = result,sig_show ="abc",ncol = 1)
  p1_2 = res[[1]]
  p1_2
  res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = c(3),result = result,sig_show ="abc",ncol = 1)
  p1_3 = res[[1]]
  p1_3

  return(list(p1_1,p1_2,p1_3,distance.data = data))
}


#' @title Calculate and Visualize Bray-Curtis Distances Between Groups
#'
#' @description
#' The `distance.micro2` function computes Bray-Curtis distances between groups of microbial community samples.
#' It generates pairwise group comparisons and visualizes the results using boxplots and bar plots with statistical annotations.
#'
#' @param ps A `phyloseq` object containing microbial community data.
#' @param group A character string specifying the grouping variable in the metadata (default: "Group").
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{p1_1}: A ggplot object showing boxplots of Bray-Curtis distances.
#'   \item \code{p1_2}: A ggplot object showing bar plots of Bray-Curtis distances.
#'   \item \code{distance.data}: A data frame with pairwise distances and corresponding group comparisons.
#' }
#'
#' @details
#' This function calculates the Bray-Curtis distances between all samples in the dataset using the `vegan` package.
#' For each pair of groups specified by the `group` variable, it extracts the distances between samples
#' and organizes the results into a data frame. The function then visualizes the distances using boxplots
#' and bar plots with statistical annotations added through the `EasyStat` package.
#'
#' @examples
#' # Example usage:
#' library(phyloseq)
#' library(EasyStat)
#' library(vegan)
#'
#' # Assuming `ps` is a phyloseq object with grouping variable "Treatment"
#' result <- distance.micro2(ps = ps, group = "Treatment")
#'
#' # Access and display the plots
#' result[[1]]  # Boxplot
#' result[[2]]  # Bar plot
#' result[[3]]  # Combined plot
#'
#' # Access the distance data
#' head(result$distance.data)
#'
#' @export

distance.micro2 <- function(ps = NULL, group = "Group") {
  map = as.data.frame(sample_data(ps))
  gro = map[, group] %>% unique()
  colnames(gro) = "group"
  conbgroup = combn(gro$group, 2)

  bray_curtis = vegan::vegdist(vegan_otu(ps), method = "bray")
  bray_curtis = as.matrix(bray_curtis)

  for (i in 1:ncol(conbgroup)) {
    a = conbgroup[, i]
    map = as.data.frame(sample_data(ps))

    chose1 = map[as.matrix(map[, group]) %>% as.vector() == a[1], ] %>% row.names()
    chose2 = map[as.matrix(map[, group]) %>% as.vector() == a[2], ] %>% row.names()

    dat = data.frame(
      group    = paste(a[1], a[2], sep = "_VS_"),
      Distance = bray_curtis[chose1, chose2] %>% as.dist() %>% as.vector()
    )

    if (i == 1) table = dat else table = rbind(table, dat)
  }

  table$id = 1:nrow(table)
  data <- table %>% dplyr::select(id, dplyr::everything())


  if (length(unique(data$group)) < 2) {
    p1_1 <- ggplot(data, aes(x = group, y = Distance)) +
      geom_boxplot() +
      theme_bw()

    p1_2 <- p1_1
    p1_3 <- p1_1

    return(list(p1_1, p1_2, p1_3, distance.data = data))
  }

  # 否则走原来的 EasyStat 流程
  result  <- EasyStat::MuiKwWlx(data = data, num = 3)
  result1 <- EasyStat::FacetMuiPlotresultBox(data = data, num = 3,
                                             result = result, sig_show = "abc", ncol = 1)
  p1_1 <- result1[[1]]

  res2 <- EasyStat::FacetMuiPlotresultBar(data = data, num = 3,
                                          result = result, sig_show = "abc", ncol = 1)
  p1_2 <- res2[[1]]

  res3 <- EasyStat::FacetMuiPlotReBoxBar(data = data, num = 3,
                                         result = result, sig_show = "abc", ncol = 1)
  p1_3 <- res3[[1]]

  return(list(p1_1, p1_2, p1_3, distance.data = data))
}




