#' @title Perform Partial Least Squares Discriminant Analysis (PLS-DA) on metabolite data
#' @description
#' The plsda.ms function conducts Partial Least Squares Discriminant Analysis (PLS-DA) on metabolite data using the mixOmics package.
#' It visualizes the results of PLS-DA including score plots, explained variance, and ellipses.

#' @param map A data frame containing sample metadata.
#' @param ps A phyloseq format file used as an alternative for the input containing metabolite composition table,
#' metabolite classification table, and sample metadata.
#' @param count A matrix of metabolite counts.
#' @param Group Column name for groupID in map table(sample metadata).
#' @return A list containing the following components:
#' \item{p}{Plot showing the PLS-DA analysis results.}
#' \item{plotdata}{Data frame containing the coordinates of samples in the PLS-DA space.}
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' library(mixOmics)
#' library(ggplot2)
#' res = plsda_ms(ps=ps.ms,Group = "Group")
#' p11 = res [[1]]
#' p11
#' dat = res[[2]]
#' dat
#' @export

plsda_ms = function(map = NULL,ps = NULL,count =NULL,
                    Group = "Group"){
count = vegan_otu(ps)
map = as.data.frame(sample_data(ps))
map$Group=factor(map$Group)
plsda.datatm <-mixOmics::plsda(count, map$Group, ncomp = 2)
#PLS-DA without centroids
mi=c("#1B9E77" ,"#D95F02")
plotIndiv(plsda.datatm , comp = c(1,2),
          group = map$Group, style = 'ggplot2' )
plotIndiv(plsda.datatm , comp = c(1,2),
          group = map$Group, style = 'ggplot2',ellipse = TRUE,
          size.xlabel = 20, size.ylabel = 20, size.axis = 25, pch = 15, cex = 5)
#----提取数据作图
a = unclass(plsda.datatm)
#--提取坐标值
plotdata = as.data.frame(a$variates$X)
plotdata$SampleType = map$Group
#-提取解释度
eig = a$prop_expl_var$X
eig[1]

# library(BiocManager)
# install("ggalt")
p = ggplot(data = plotdata,aes(x=comp1,y=comp2,group=SampleType,color=SampleType))+geom_point(size=5)+
  stat_ellipse(type = "t", linetype = 2)+
  geom_encircle(s_shape=1, expand=0) +
  labs(x=paste("X-variate 1 (", format(round(100 * eig[1],2)), "%)", sep=""),
       y=paste("X-variate 2 (", format(round(100 * eig[2],2)), "%)", sep=""))+
  labs(title = "PLS-DA")
p

mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A")
p=p+theme_bw()+scale_colour_manual(values = mi)+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  geom_hline(yintercept=0) + geom_vline(xintercept=0)
p
return(list(p,plotdata))
}



