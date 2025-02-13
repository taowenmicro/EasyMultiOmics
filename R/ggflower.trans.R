#' @title Create a circular flower plot based on transcriptome functional composition data
#' @description
#'  The ggflower.metf function generates a circular flower plot based on transcriptome functional composition data.
#'  The plot visualizes the presence or absence of genes
#'  across different sample groups in a circular layout.
#' @param otu Transcriptome functional composition table.
#' @param tax Transcriptome functional classification table.
#' @param map Sample metadata.
#' @param ps A phyloseq format file used as an alternative for the input containing transcriptome functional composition table, tax, and sample metadata.
#' @param group Column name for groupID in map table(sample metadata).
#' @param rep Repeat number of transcriptome functional composition data.
#' @param m1 Petal shape, square to round to prismatic, the value gradually decreases.
#' @param start The rotation angle of the petals, the greater the value, the greater the angle.
#' @param a The width of the petals.
#' @param b Distance from petal to center.
#' @param lab.leaf The distance from the label to the center of the circle.
#' @param col.cir Center color.
#' @return A list containing the following components:
#' \item{p}{Circular flower plot visualizing the presence or absence of genes across sample groups.}
#' \item{ven2}{Processed genes data frame used for plotting.}
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' library(phyloseq)
#' library(ggplot2)
#' data(ps.trans)
#' ps = ps.trans%>% filter_taxa(function(x) sum(x ) > 5 , TRUE)
#' res <- ggflower.trans(ps= ps,group = "Group",start = 1,
#' m1 = 1,a = 0.2, b = 1,lab.leaf = 1col.cir = "yellow",N = 0.5 )
#' p14 = res[[1]]
#' p14
#' dat = res[[2]]
#' dat
#' @export


ggflower.trans = function(otu = NULL,tax = NULL,map = NULL,ps = NULL,
                         group = "Group",
                         rep = 6,
                         m1 = 2,
                         start = 1,
                         a = 0.2,
                         b = 1,
                         lab.leaf = 1,
                         col.cir = "yellow",
                         a.cir = 0.5,
                         b.cir = 0.5,
                         m1.cir = 2,
                         N = 0.5
){

  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)

  mapping = as.data.frame(phyloseq::sample_data(ps))

  aa = ggClusterNet::vegan_otu(ps)
  otu_table = as.data.frame(t(aa))
  count = aa
  sub_design <- as.data.frame(phyloseq::sample_data(ps))
  sub_design$SampleType = sub_design$Group
  phyloseq::sample_data(ps ) = sub_design
  count[count > 0] <- 1
  count2 = as.data.frame(count)
  # group
  iris.split <- split(count2,as.factor(sub_design$Group))
  #group mean
  iris.apply <- lapply(iris.split,function(x)colSums(x[]))
  # conbine result
  iris.combine <- do.call(rbind,iris.apply)
  ven2 = t(iris.combine)

  for (i in 1:length(unique(phyloseq::sample_data(ps)$Group))) {
    aa <- as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,1]
    bb =  as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,2]
    ven2[,aa] = ven2[,aa]/bb
  }

  ven2[ven2 < N]  = 0
  ven2[ven2 >=N]  = 1
  ven2 = as.data.frame(ven2)
  ven3 = as.list(ven2)
  ven2 = as.data.frame(ven2)
  all_num = dim(ven2[rowSums(ven2) == length(levels(sub_design$Group)),])[1]
  ven2[,1] == 1
  A = rep("A",length(colnames(ven2)))
  B = rep(1,length(colnames(ven2)))
  i = 1
  for (i in 1:length(colnames(ven2))) {
    B[i] = length(ven2[rowSums(ven2) == 1,][,i][ven2[rowSums(ven2) == 1,][,i] == 1])
    A[i] = colnames(ven2)[i]
  }
  n   <- length(A)
  deg <- 360 / n
  t = 1:n
  print(deg)

  p <- ggplot() +
    # geom_point(aes(x = 5 + cos((start + deg * (t - 1)) * pi / 180) * lab.leaf, y = 5 + sin((start + deg * (t - 1)) * pi / 180) *lab.leaf)) +
    ggforce::geom_ellipse(aes(x0 = 5 + cos((start + deg * (t - 1)) * pi / 180),
                              y0 = 5 + sin((start + deg * (t - 1)) * pi / 180),
                              a = a,
                              b = b,
                              angle = (n/2 +seq(0,1000,2)[1:n])/n * pi,
                              m1 = m1,
                              fill = as.factor(1:n)),show.legend = F) +
    ggforce::geom_ellipse(aes(x0 = 5,y0 = 5,a = a.cir,b = b.cir,angle = 0,m1 = m1.cir),fill = col.cir) +
    geom_text(aes(x = 5,y = 5,label = paste("OVER :",all_num,sep = ""))) +
    geom_text(aes(
      x = 5 + cos((start + deg * (t - 1)) * pi / 180) * lab.leaf,
      y = 5 + sin((start + deg * (t - 1)) * pi / 180) * lab.leaf,
      label = paste(A,":",B,sep = "")),angle = 360/n*((1:n)-1)  ) +
    coord_fixed() + theme_void()
  p
  return(list(p ,ven2))
}




