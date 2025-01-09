#' @title Generate Circular Dendrogram with Bar Plots for Microbial Data
#'
#' @description
#' This function generates a circular hierarchical clustering dendrogram combined with bar plots of microbial composition.
#' It uses taxonomic relative abundances at the Phylum level (default) and allows grouping of samples into clusters.
#'
#' @param ps A `phyloseq` object containing microbiome data.
#' @param Top An integer specifying the number of top taxa to display. Taxa not in the top are grouped as `"Other"`. Default is `15`.
#' @param dist A character string specifying the distance metric for clustering. Default is `"bray"`.
#' @param cuttree An integer specifying the number of clusters to cut the dendrogram into. Default is `3`.
#' @param hcluter_method A character string specifying the hierarchical clustering method. Default is `"complete"`.
#'
#' @return
#' A list containing the following elements:
#' \describe{
#'   \item{A `ggtree` circular dendrogram with bar plots showing relative abundances of taxa.}
#'   \item{A data frame containing the processed data used for plotting, including taxonomic annotations and abundances.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Standardizes the input microbiome data using `ggClusterNet::scale_micro`.
#'   \item Calculates a distance matrix using the specified method (`dist`) and performs hierarchical clustering.
#'   \item Cuts the clustering dendrogram into the specified number of clusters (`cuttree`) and combines cluster information with sample metadata.
#'   \item Groups microbial taxa at the Phylum level and calculates relative abundances.
#'   \item Selects the top `Top` taxa by abundance, grouping the rest as `"Other"`.
#'   \item Combines the circular dendrogram with bar plots showing the relative abundances of taxa for each sample.
#' }
#'
#' @examples
#' \dontrun{
#' res = cir_barplot.metm(ps = ps.16s,Top = 15,dist = "bray",cuttree = 3,hcluter_method = "complete")
#' p17 = res[[1]]
#' dat= res[[2]]
#'}
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#'
#' @export
cir_barplot.metm = function(
    ps = ps,
    Top = 15,
    dist = "bray",
    cuttree = 3,
    hcluter_method = "complete"
){
  # phyloseq(ps)对象标准化
  otu = ps %>%
    ggClusterNet::scale_micro() %>%
    ggClusterNet::vegan_otu() %>% t() %>%
    as.data.frame()
  # 导出OTU表
  # otu = as.data.frame(t(vegan_otu(ps1_rela)))

  #计算距离矩阵
  unif = phyloseq::distance(ps %>% ggClusterNet::scale_micro() , method = dist)
  # 聚类树，method默认为complete
  hc <- stats::hclust(unif, method = hcluter_method )
  #  take grouping with hcluster tree
  clus <- cutree(hc, cuttree )

  d = data.frame(label = names(clus),
                 member = factor(clus))
  # eatract mapping file
  map = as.data.frame(phyloseq::sample_data(ps))

  dd = merge(d,map,by = "row.names",all = F)
  row.names(dd) = dd$Row.names
  dd$Row.names = NULL
  # library(tidyverse)
  # ggtree绘图 #----
  p  = ggtree::ggtree(hc, layout='circular') %<+% dd +
    geom_tippoint(size=5, shape=21, aes(fill= Group, x=x)) +
    geom_tiplab(aes(color = Group,x=x * 1.2), hjust=1,offset=0.3) + xlim(-0.5,NA)
  p

  psdata =  ggClusterNet::tax_glom_wt(ps = ps %>% ggClusterNet::scale_micro(),ranks = "Phylum" )
  # 转化丰度值
  psdata = psdata%>% phyloseq::transform_sample_counts(function(x) {x/sum(x)} )

  #--提取otu和物种注释表格
  otu = phyloseq::otu_table(psdata)
  tax = phyloseq::tax_table(psdata)
  head(tax)
  #--按照指定的Top数量进行筛选与合并
  j = "Phylum"
  for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {
      tax[i,j] =tax[i,j]
    } else {
      tax[i,j]= "Other"
    }
  }
  phyloseq::tax_table(psdata)= tax

  Taxonomies <- psdata %>% phyloseq::psmelt()

  Taxonomies$Abundance = Taxonomies$Abundance * 100

  Taxonomies$OTU = NULL
  colnames(Taxonomies)[1] = "id"

  head(Taxonomies)

  dat2 = data.frame(id = Taxonomies$id,Abundance = Taxonomies$Abundance,Phylum = Taxonomies$Phylum)
  head(dat2)

  p2 <- p +
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(
      data=dat2,
      geom=geom_bar,
      mapping=aes(
        x=Abundance,
        y=id,
        fill= Phylum
      ),
      stat="identity",
      width = 0.4,
      orientation="y",
      offset=0.9,
      pwidth=2,
      axis.params=list(
        axis = "x",
        text.angle = -45,
        hjust = 0,
        vjust = 0.5,
        nbreak = 4
      )
    )  +theme_void()

  return(list(plot = p2,plotdata = dat2))

}

