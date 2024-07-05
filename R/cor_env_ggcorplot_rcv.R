cor_env_ggcorplot_rcv = function(ps =  ps.tem,
                                 jj = NULL,
                                 tran  =TRUE,
                                 Top = 10,
                                 env1 = ftab,
                                 label =  TRUE,
                                 col_cluster = TRUE,
                                 row_cluster = TRUE,
                                 method = "spearman",
                                 r.threshold= 0,
                                 p.threshold= 0,
                                 theme.size = 10){

  if (is.null(jj)) {
    ps_tem = ps %>%
      ggClusterNet::scale_micro(method = "rela")
  } else{
    ps_tem = ps%>%
      ggClusterNet::scale_micro(method = "rela") %>%
      ggClusterNet::tax_glom_wt(ranks = j)
  }



  id <- ps_tem %>%
    ggClusterNet::scale_micro(method = "rela") %>%
    # ggClusterNet::filter_OTU_ps(100) %>%
    ggClusterNet::vegan_otu() %>%
    t() %>% as.data.frame() %>%
    rowCV %>%
    sort(decreasing = TRUE) %>%
    head(Top) %>%
    names()

  otu = phyloseq::otu_table(ps_tem)
  tax = phyloseq::tax_table(ps_tem)

  head(otu)

  data = otu[id,] %>% t() %>%as.data.frame()

  env2 = data
  env2 <- env2[match(row.names(env1),row.names(env2)),]
  env0 <- cbind(env1,env2)
  occor = psych::corr.test(env0,use="pairwise",method=method,adjust="fdr",alpha=.05)
  occor.r = occor$r
  occor.p = occor$p

  occor.r[occor.p > p.threshold&abs(occor.r) < r.threshold] = 0

  head(env0)
  # data[data > 0.3]<-0.3
  #drop gene column as now in rows

  if (col_cluster ==  TRUE) {
    clust <- hclust(dist(env1 %>% as.matrix()%>% t())) # hclust with distance matrix
    ggtree_plot <- ggtree::ggtree(clust)
  }
  if (row_cluster ==  TRUE) {
    v_clust <- hclust(dist(env2 %>% as.matrix() %>% t()))
    ggtree_plot_col <- ggtree::ggtree(v_clust) + ggtree::layout_dendrogram()
  }

  head(occor.r)
  occor.r = as.data.frame(occor.r)
  data <- occor.r[colnames(env1),colnames(env2)]
  data$id = row.names(data)
  pcm = reshape2::melt(data, id = c("id"))

  head(pcm)

  p1 = ggplot(pcm, aes(y = id, x = variable)) +
    # geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) +
    geom_tile(aes(size = value,fill = value))+
    scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) +
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  +
    # scale_fill_manual(values = colours, guide = FALSE) +
    scale_x_discrete(limits = rev(levels(pcm$variable)))  +
    scale_y_discrete(position = "right") +
    scale_fill_gradientn(colours =colorRampPalette(c("#377EB8","#F7F4F9","#E41A1C"))(60)) +
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.x = element_text(colour = "black",size = theme.size,angle = 60,vjust = 1,hjust = 1)
    )

  # colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
  #              "#F6AE2D","#86BBD8")
  #----样本在y轴上
  p2 = ggplot(pcm, aes(y = id, x = variable)) +
    geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) +
    scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) +
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  +
    # scale_fill_manual(values = colours, guide = FALSE) +
    scale_x_discrete(limits = rev(levels(pcm$variable)))  +
    scale_y_discrete(position = "right")  +
    scale_fill_gradientn(colours =colorRampPalette(c("#377EB8","#F7F4F9","#E41A1C"))(60))  +
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.x = element_text(colour = "black",size = theme.size,angle = 60,vjust = 1,hjust = 1)
    )

  if (col_cluster ==  TRUE) {
    p1 <- p1  %>%
      aplot::insert_left(ggtree_plot, width=.2)
    p2 <- p2  %>%
      aplot::insert_left(ggtree_plot, width=.2)
  }

  if (label ==  TRUE) {
    p1 <- p1  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
    p2 <- p2  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
  }

  return(list(p1,p2,env2,data))


}
