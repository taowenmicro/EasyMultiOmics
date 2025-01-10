
#' @title FEAST Source Contribution Visualization
#'
#' @description
#' The `Plot_FEAST` function generates a polar bar chart (pie chart) to visualize the relative contributions of different sources using FEAST (Fast Expectation-Maximization for Source Tracking).
#'
#' @param data A matrix or data frame where rows represent sources, and columns represent contributions across samples. Typically, this is the output from FEAST analysis. Default is `result`.
#'
#' @return A ggplot2 object representing the polar bar chart with source contributions and their respective percentages.
#'
#' @details
#' The function calculates the mean contribution of each source across samples, normalizes the contributions to ensure they sum to 100%, and visualizes the contributions as a polar bar chart (pie chart). Each slice of the chart represents a source, and the label includes the source name and its percentage contribution.
#'
#' @examples
#' \dontrun{
#' result = FEAST.micro(ps = ps.16s,group = "Group",sinkG = "OE",sourceG = c("WT","KO"))
#' p <- Plot_FEAST(data = result)
#' }
#'
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#'
#' @export


Plot_FEAST = function(data = result){

  asx = as.data.frame(rowMeans(result))

  asx  = as.matrix(asx)
  asx_norm = t(t(asx)/colSums(asx)) #* 100 # normalization to total 100
  head(asx_norm)
  asx_norm = as.data.frame( asx_norm)
  colnames(asx_norm) = "present"
  # plotname = paste(path,"/FEAST_mean.pdf",sep = "")
  # pdf(file = plotname,width = 6,height = 6)
  labs <- paste0(row.names(asx_norm)," (", round(asx_norm[,1]/sum(asx_norm[,1])*100,2), "%)")
  asx_norm$lab = labs
  asx_norm$ID = row.names( asx_norm)
  head(asx_norm)
  p <-  ggplot(asx_norm, aes( x = "",y = present, fill = labs)) +
    geom_bar(stat = "identity",width = 20,color = "black") +
    coord_polar(theta = "y",direction=1) +
    theme_void()
  return(p)
}
