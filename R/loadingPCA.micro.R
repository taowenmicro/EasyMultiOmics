#' @title PCA loading Matrix screening of characteristic microorganisms
#' @description
#' This function conducts PCA analysis on microbiome data to
#' extract the loading Matrix, screen for characteristic microorganisms and visualize them.
#' The importance of the variables is  sorted according to
#' the square value of the correlation between the variable and the PC1 axis.
#' @param ps A phyloseq format file used as an alternative for the input containing otu, tax, and map.
#' @param Top The top microorganisms to visualize.
#' @returns A list object containing the following components:
#' \item{p}{A PCA correlation plot of the selected number of characteristic microorganisms and
#'  the correlation decreases from top to bottom.}
#' \item{index}{Data frame containing the PCA load matrix and relative abundance of all microorganisms.}
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' res = loadingPCA.micro(ps = ps.16s,Top = 20)
#' p34.1 = res[[1]]
#' p34.1
#' dat = res[[2]]
#' dat
loadingPCA.micro = function(ps = ps,Top = 20){
  count =ps %>%
    ggClusterNet::vegan_otu() %>% t()
  count[is.na(count)] = 0
  norm = t(t(count)/colSums(count,na=TRUE))# * 100 # normalization to total 100
  otu.pca <- stats::prcomp(t(norm), scale= FALSE)

  #提取
  yangpin<-otu.pca$x
  yangpin=as.data.frame(yangpin)
  yangpin$SampleType= phyloseq::sample_data(ps)$Group
  #提取荷载坐标
  bianliang<-otu.pca$rotation
  bianliang=as.data.frame(bianliang)
  head(bianliang)
  dim(norm)
  index = merge(norm ,bianliang, by="row.names",all=F)
  head(index)
  row.names(index)=index$Row.names
  index$Row.names=NULL
  head(index)
  index$id = row.names(index)
  ##手动选择10个最终要的变量 PCA载荷矩阵挑选37个成分提取差异.txt
  index$PCone = index$PC1^2
  top = index %>% arrange(desc(PCone)) %>%head(Top)
  head(top)
  top$ID  = top$id
  head(top)

  p=ggplot(top, aes(x = PCone, y = reorder(ID,PCone)))  +
    geom_segment(aes(yend=ID),xend=0,size=3,colour = "#1B9E77" )+
    geom_point(size=4,pch=20, colour = "#1B9E77")+theme_bw()+
    theme(axis.text.x = element_text(colour = "black",size = 20,face = "bold"),
          axis.text.y = element_text(colour = "black",size = 10,face = "bold"))
  p
  return(list(p,index))
}


