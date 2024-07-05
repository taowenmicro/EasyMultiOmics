


# res = cluster_plot (ps= ps,hcluter_method = "complete",
#   dist = "bray",cuttree = gnum,row_cluster = T,col_cluster =  T)

cluster.omics = function(
    ps= ps,
    hcluter_method = "complete",
    dist = "bray",
    cuttree = gnum,
    row_cluster = T,
    col_cluster =  T

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
  dd = merge(d,map,by = "row.names",all = FALSE)
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


