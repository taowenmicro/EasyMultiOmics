#' @title Cluster and Visualize transcriptome functional composition
#' @description
#' This function performs hierarchical clustering and visualizes transcriptome functional composition
#' using phyloseq data. It provides various visualizations, including clustering dendrograms
#' and heatmaps, to explore relationships among samples.
#'
#' @param ps A `phyloseq` object containing transcriptome functional composition data.
#' @param hcluter_method A string specifying the hierarchical clustering method.
#' Default is `"complete"`. Options include `"single"`, `"average"`, and `"ward.D"`.
#' @param dist A string specifying the distance method for clustering.
#' Default is `"bray"`. Other options include methods available in `phyloseq::distance`.
#' @param cuttree An integer specifying the number of clusters for cutting the dendrogram.
#' @param row_cluster Logical. Whether to perform clustering on rows. Default is `TRUE`.
#' @param col_cluster Logical. Whether to perform clustering on columns. Default is `TRUE`.
#' @return A list containing the following elements:
#' \item{p0}{A ggtree plot of the clustering dendrogram.}
#' \item{p1}{A heatmap showing relative abundances with optional row and column clustering.}
#' \item{p2}{A heatmap with sample IDs on the y-axis, showing relative abundances.}
#' \item{tem}{A data frame representing the distance matrix used for clustering.}
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' data(ps.trans)
#' ps = ps.trans%>% filter_taxa(function(x) sum(x ) > 5 , TRUE)
#' res = cluster.trans(ps= ps,hcluter_method = "complete",dist = "bray",cuttree = 3,row_cluster = TRUE,col_cluster =  TRUE)
#' p4 = res[[1]]
#' p4
#' p4_1 = res[[2]]
#' p4_1
#' p4_2 = res[[3]]
#' p4_2
#' dat = res[[4]]# cluster distance
#' head(dat)
cluster_trans = function(
    ps= ps,
    hcluter_method = "complete",
    dist = "bray",
    cuttree = gnum,
    row_cluster = TRUE,
    col_cluster =  TRUE

){
  # phyloseq(ps)对象标准化
  ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) )
  # 导出OTU表
  otu = as.data.frame(t(ggClusterNet::vegan_otu(ps1_rela)))

  #计算距离矩阵
  unif = phyloseq::distance(ps1_rela , method = dist)
  # 聚类树，method默认为complete
  hc <- stats::hclust(unif, method = hcluter_method )

  #  take grouping with hcluster tree
  clus <- stats::cutree(hc, cuttree )
  # 提取树中分组的标签和分组编号
  d = data.frame(label = names(clus),
                 member = factor(clus))
  # eatract mapping file
  map = as.data.frame(phyloseq::sample_data(ps))
  # 合并树信息到样本元数据
  dd = merge(d,map,by = "row.names",all = F)
  row.names(dd) = dd$Row.names
  dd$Row.names = NULL

  library(ggtree)
  # ggtree绘图#----
  p0  = ggtree::ggtree(hc) %<+% dd +
    geom_tippoint(size=5, shape=21, aes(fill= Group, x=x)) +
    geom_tiplab(aes(color = Group,x=x * 1.2), hjust=1)
  p0


  if (col_cluster ==  TRUE) {
    clust <- hclust(unif) # hclust with distance matrix
    ggtree_plot <- ggtree::ggtree(clust)
  }
  if (row_cluster ==  TRUE) {
    clust <- hclust(unif)
    ggtree_plot_col <- ggtree::ggtree(clust) + ggtree::layout_dendrogram()
  }


  tem = unif %>% as.matrix()
  tem = 1- tem
  tem = tem %>% as.data.frame()
  tem$id = row.names(tem)
  pcm = reshape2::melt(tem, id = c("id"))
  head(pcm)
  # pcm$variable = factor(pcm$variable,levels = map$ID)
  # pcm$id = factor(pcm$id,levels = rig$id)

  p1 = ggplot(pcm, aes(y = id, x = variable)) +
    # geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) +
    geom_tile(aes(fill = value))+
    scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) +
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  +
    # scale_fill_manual(values = colours, guide = FALSE) +
    scale_x_discrete(limits = (levels(pcm$variable)))  +
    scale_y_discrete(limits = rev(levels(pcm$variable)),position = "right") +
    scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral")[11:1])(60))+
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.x = element_text(colour = "black",angle = 90,hjust = 1,
                                 vjust = 0)

    )
  p1
  #----样本在y轴上
  p2 = ggplot(pcm, aes(y = id, x = variable)) +
    geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) +
    scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) +
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  +
    # scale_fill_manual(values = colours, guide = FALSE) +
    scale_x_discrete(limits = (levels(pcm$variable)))  +
    scale_y_discrete(limits = rev(levels(pcm$variable)),position = "right") +
    scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral")[11:1])(60)) +
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.x = element_text(colour = "black",angle = 90,hjust = 1,
                                 vjust = 0)

    )
  p2



  if (col_cluster ==  TRUE) {
    p1 <- p1  %>%
      aplot::insert_left(ggtree_plot, width=.2)
    p2 <- p2  %>%
      aplot::insert_left(ggtree_plot, width=.2)
  }



  if (row_cluster ==  TRUE) {
    p1 <- p1  %>%
      aplot::insert_top(ggtree_plot_col, height=.2)
    p2 <- p2  %>%
      aplot::insert_top(ggtree_plot_col, height=.2)
  }
  return(list(p0,p1,p2,tem))
}


