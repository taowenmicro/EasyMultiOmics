
heatmap_GSVA <- function(heat = sub_heat,
                         map = sub_map,
                         col_cluster =  F,
                         row_cluster =  F,
                         label =  TRUE){
  
  heat = as.data.frame(heat)
  heat$id = row.names(heat)
  head(heat)
  data = heat %>% 
    dplyr::select(id,everything())
  # data[data > 0.3]<-0.3
  mat <- data[,-1] #drop gene column as now in rows
  
  if (col_cluster) {
    clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix
    ggtree_plot <- ggtree::ggtree(clust)
  }
  if (row_cluster) {
    v_clust <- hclust(dist(mat %>% as.matrix() %>% t()))
    ggtree_plot_col <- ggtree::ggtree(v_clust) + ggtree::layout_dendrogram()
  }
  map = map %>%
    as.tibble() %>%
    arrange(desc(Group))
  map$ID = factor(map$ID,levels = unique(map$ID))
  if (label ) {
    
    
    labels = ggplot(map, aes(x = ID, y=1, fill=Group)) + geom_tile() +
      scale_fill_brewer(palette = 'Set1',name="Cell Type") +
      theme_void()
  }
  
  pcm = reshape2::melt(data, id = c("id"))
  head(pcm)
  
  pcm$variable = factor(pcm$variable,levels = as.character(map$ID))
  
  
  p1 = ggplot(pcm, aes(y = id, x = variable)) + 
    # geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) + 
    geom_tile(aes(size = value,fill = value))+
    scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) + 
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  + 
    # scale_fill_manual(values = colours, guide = FALSE) + 
    scale_x_discrete(limits = rev(levels(pcm$variable)))  + 
    scale_y_discrete(position = "right") +
    scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral")[11:1])(60))+
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.x = element_text(colour = "black",angle = 90)
      
    )
  
  colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
               "#F6AE2D","#86BBD8")
  #----样本在y轴上
  p2 = ggplot(pcm, aes(y = id, x = variable)) + 
    geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) + 
    scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) + 
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  + 
    # scale_fill_manual(values = colours, guide = FALSE) + 
    scale_x_discrete(limits = rev(levels(pcm$variable)))  + 
    scale_y_discrete(position = "right")  +
    scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral")[11:1])(60)) +
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.x = element_text(colour = "black",angle = 90)
      
    )
  
  if (col_cluster) {
    p1 <- p1  %>%
      aplot::insert_left(ggtree_plot, width=.2) 
    p2 <- p2  %>%
      aplot::insert_left(ggtree_plot, width=.2) 
  }
  
  if (label) {
    p1 <- p1  %>%
      aplot::insert_top(labels, height=.02) 
    p2 <- p2  %>%
      aplot::insert_top(labels, height=.02) 
  }
  
  if (row_cluster) {
    p1 <- p1  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
    p2 <- p2  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
  }
  
  return(list(p1,p2))
  
  
}
