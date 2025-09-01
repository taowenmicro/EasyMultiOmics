#' @title Multiple FEAST Source Contribution Visualizations
#'
#' @description
#' The `MuiPlot_FEAST` function creates multiple polar bar charts (pie charts) to visualize the source contributions for each sample or group in the result matrix from FEAST analysis. Each chart represents the relative contributions of sources to a particular sample/group.
#'
#' @param data A matrix or data frame where rows represent sources, and columns represent samples or groups. Typically, this is the output from FEAST analysis. Default is `result`.
#'
#' @return A combined ggplot2 object with polar bar charts for each sample/group. Each chart includes source contributions and their respective percentages.
#'
#' @details
#' The function iterates through the columns of the input matrix, calculates the mean contribution of each source for a given sample/group, normalizes the contributions, and visualizes them as polar bar charts. The charts are combined into a single plot with a shared legend for easy comparison.
#'
#' @examples
#' \dontrun{
#' result = FEAST.micro(ps = ps.16s,group = "Group",sinkG = "OE",sourceG = c("WT","KO"))
#' p2 = MuiPlot_FEAST(data = result)
#' }
#'
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#'
#' @export

MuiPlot_FEAST = function(data = result){

  par(mfrow=c(2,dim(result)[2]/2), mar=c(1,1,1,1))
  # layouts = as.character(unique(metadata$SampleType))
  i = 1
  plots = list()
  for (i in 1:length(colnames(result))) {

    asx = data.frame(row.names = row.names(result),result[,i])

    asx  = as.matrix(asx)
    asx_norm = t(t(asx)/colSums(asx)) #* 100 # normalization to total 100
    head(asx_norm)
    asx_norm = as.data.frame( asx_norm)
    colnames(asx_norm) = colnames(result)[i]
    labs <- paste0(row.names(asx_norm)," (", round(asx_norm[,1]/sum(asx_norm[,1])*100,2), "%)")
    asx_norm$lab = labs
    asx_norm$ID = row.names( asx_norm)


    x <- colnames(result)[i]
    p <-  ggplot(asx_norm, aes( x = "",y = !!sym(x), fill = lab)) +
      geom_bar(stat = "identity",width = 20,color = "black") +
      coord_polar(theta = "y",direction=1) +
      labs(title = x) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(fill = guide_legend(title = NULL))
    p
    plots[[i]] = p
  }

  p  = ggpubr::ggarrange(plotlist = plots, common.legend = TRUE, legend="right")
  return(p)
}

