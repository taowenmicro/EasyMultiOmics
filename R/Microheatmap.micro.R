



Microheatmap.micro <- function(ps_rela,
                         id,
                         label =  TRUE,
                         col_cluster = TRUE,
                         row_cluster = TRUE,
                         ord.col = TRUE,
                         scale = TRUE,# 是否标准化丰度，可以避免极大丰度的影响
                         axis_order.s = axis_order.s,
                         row.lab = NULL,
                         col1 = ggsci::pal_gsea(alpha = 1)(12)
                         ){

  map = phyloseq::sample_data(ps_rela)
  map$ID = row.names(map)
  phyloseq::sample_data(ps_rela) = map
  otu = as.data.frame(t(ggClusterNet::vegan_otu(ps_rela)))
  otu = as.matrix(otu[id,])
  ps_heatm = phyloseq::phyloseq(
    phyloseq::otu_table(otu,taxa_are_rows = TRUE),
    phyloseq::tax_table(ps_rela),
    phyloseq::sample_data(ps_rela)

  )

  # print(ps_heatm)
  datah <- as.data.frame(t(ggClusterNet::vegan_otu(ps_heatm)))
  head(datah)
  tax = as.data.frame(ggClusterNet::vegan_tax(ps_heatm))

  otutaxh = cbind(datah,tax)
  head(otutaxh)

  otutaxh$id = paste(row.names(otutaxh))
  # otutaxh$id  = row.names(otutaxh)
  row.names(otutaxh) = otutaxh$id


  data <- otutaxh[,c("id",phyloseq::sample_names(ps_rela))]

  rig <- data[,phyloseq::sample_names(ps_rela)] %>% rowMeans() %>% as.data.frame()
  head(rig)
  colnames(rig) = "MeanAbundance"
  rig$id = row.names(rig)
  rig = rig %>% dplyr::arrange(MeanAbundance)
  rig$id = factor(rig$id,levels = rig$id)
  head(rig)
  p_rig = ggplot(rig) + geom_bar(aes(y = id,x = MeanAbundance),
                                 fill = "#A54657",
                                 stat = "identity") + theme_void()

  tem = data[,phyloseq::sample_names(ps_rela)] %>% as.matrix()

  if (scale == TRUE) {
    tem = scale(t(tem)) %>% t() %>%
      as.data.frame()
  } else if (scale == FALSE){
    tem = tem
  }

   data[,phyloseq::sample_names(ps_rela)] = tem


  # data[data > 0.3]<-0.3
  mat <- data[,-1] #drop gene column as now in rows

  if (row_cluster ==  TRUE) {
    clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix
    ggtree_plot <- ggtree::ggtree(clust)
  }
  if (col_cluster ==  TRUE) {
    v_clust <- hclust(dist(mat %>% as.matrix() %>% t()))
    ggtree_plot_col <- ggtree::ggtree(v_clust) + ggtree::layout_dendrogram()
  }

  if (label ==  TRUE) {
    map = as.data.frame(phyloseq::sample_data(ps_rela))
    map$ID = row.names(map)
    labels= ggplot(map, aes(x = ID, y=1, fill=Group)) + geom_tile() +
      scale_fill_brewer(palette = 'Set1',name="Cell Type") +
      theme_void()
  }

  if (!is.null(row.lab) ) {
    tax = ps_rela %>% vegan_tax() %>%
      as.data.frame()
    head(tax)
    tax$ID = row.names(tax)
    p.row.lab= ggplot(tax, aes(y = ID, x=1, fill=!!sym(row.lab))) + geom_tile() +
      scale_fill_brewer(palette = 'Set3',name="Cell Type") +
      theme_void()
  }



  map = phyloseq::sample_data(ps_rela) %>% as.tibble() %>%
    dplyr::arrange(Group) %>% as.data.frame()
  map$ID

  # tax = ps_rela %>% vegan_tax() %>% as.data.frame()
  # data$id =  tax[,new.id]
  pcm = reshape2::melt(data, id = c("id"))
  pcm$variable = factor(pcm$variable,levels = map$ID)

  if (ord.col == TRUE) {
    pcm$variable = factor(pcm$variable,levels = axis_order.s)

  }

  pcm$id = factor(pcm$id,levels = rig$id)
  head(pcm)
  # col0 =colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral")[11:1])(60)
  # col1 = ggsci::pal_gsea()(12)

  p1 = ggplot(pcm, aes(y = id, x = variable)) +
    # geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) +
    geom_tile(aes(size = value,fill = value))+
    # scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) +
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  +
    # scale_fill_manual(values = colours, guide = FALSE) +
    scale_x_discrete(limits = rev(levels(pcm$variable)))  +
    scale_y_discrete(position = "right") +
    scale_fill_gradientn(colours =col1)+
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.y = element_text(size = 3),
      axis.text.x = element_text(colour = "black",angle = 90)

    )

  colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
               "#F6AE2D","#86BBD8")
  #----样本在y轴上
head(pcm)
  p2 = ggplot(pcm, aes(y = id, x = variable)) +
    geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) +
    scale_size_continuous(limits = c(0.001, 1000), range = c(2,25), breaks = c(0.1,0.5,1)) +
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  +
    # scale_fill_manual(values = colours, guide = FALSE) +
    scale_x_discrete(limits = rev(levels(pcm$variable)))  +
    scale_y_discrete(position = "right")  +
    scale_fill_gradientn(colours = col1) +
     theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.x = element_text(colour = "black",angle = 90)

      )

  p1 <- p1  %>%
    aplot::insert_right(p_rig, width=.2)
  p2 <- p2  %>%
    aplot::insert_right(p_rig, width=.2)

  if (!is.null(row.lab)) {
    p1 <- p1  %>%
      aplot::insert_left(p.row.lab, width=.02)
    p2 <- p2  %>%
      aplot::insert_left(p.row.lab, width=.02)
  }

  if (row_cluster ==  TRUE) {
    p1 <- p1  %>%
      aplot::insert_left(ggtree_plot, width=.2)
    p2 <- p2  %>%
      aplot::insert_left(ggtree_plot, width=.2)
  }






  if (label ==  TRUE) {
    p1 <- p1  %>%
      aplot::insert_top(labels, height=.02)
    p2 <- p2  %>%
      aplot::insert_top(labels, height=.02)
  }

  if (col_cluster ==  TRUE) {
    p1 <- p1  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
    p2 <- p2  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
  }

  return(list(p1,p2,plotdata = pcm))
}


Microheatmap2 <- function(ps_rela,
                         id,
                         label =  TRUE,
                         col_cluster = TRUE,
                         row_cluster = TRUE,
                         ord.col = TRUE,
                         axis_order.s = axis_order.s,
                         row.lab = NULL
){

  map = phyloseq::sample_data(ps_rela)
  map$ID = row.names(map)
  phyloseq::sample_data(ps_rela) = map
  otu = as.data.frame(t(ggClusterNet::vegan_otu(ps_rela)))
  otu = as.matrix(otu[id,])
  ps_heatm = phyloseq::phyloseq(
    phyloseq::otu_table(otu,taxa_are_rows = TRUE),
    phyloseq::tax_table(ps_rela),
    phyloseq::sample_data(ps_rela)

  )

  # print(ps_heatm)
  datah <- as.data.frame(t(ggClusterNet::vegan_otu(ps_heatm)))
  head(datah)
  tax = as.data.frame(ggClusterNet::vegan_tax(ps_heatm))

  otutaxh = cbind(datah,tax)
  head(otutaxh)

  otutaxh$id = paste(row.names(otutaxh))
  # otutaxh$id  = row.names(otutaxh)
  row.names(otutaxh) = otutaxh$id


  data <- otutaxh[,c("id",phyloseq::sample_names(ps_rela))]

  rig <- data[,phyloseq::sample_names(ps_rela)] %>% rowMeans() %>% as.data.frame()
  head(rig)
  colnames(rig) = "MeanAbundance"
  rig$id = row.names(rig)
  rig = rig %>% dplyr::arrange(MeanAbundance)
  rig$id = factor(rig$id,levels = rig$id)
  p_rig = ggplot(rig) + geom_bar(aes(y = id,x = MeanAbundance),
                                 fill = "#A54657",
                                 stat = "identity") + theme_void()

  tem = data[,phyloseq::sample_names(ps_rela)] %>% as.matrix()

  # tem = scale(t(tem)) %>%
  #   t() %>%
  #   as.data.frame()
  data[,phyloseq::sample_names(ps_rela)] = tem


  # data[data > 0.3]<-0.3
  mat <- data[,-1] #drop gene column as now in rows

  if (row_cluster ==  TRUE) {
    clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix
    ggtree_plot <- ggtree::ggtree(clust)
  }
  if (col_cluster ==  TRUE) {
    v_clust <- hclust(dist(mat %>% as.matrix() %>% t()))
    ggtree_plot_col <- ggtree::ggtree(v_clust) + ggtree::layout_dendrogram()
  }

  if (label ==  TRUE) {
    map = as.data.frame(phyloseq::sample_data(ps_rela))
    map$ID = row.names(map)
    labels= ggplot(map, aes(x = ID, y=1, fill=Group)) + geom_tile() +
      scale_fill_brewer(palette = 'Set1',name="Cell Type") +
      theme_void()
  }

  if (!is.null(row.lab) ) {
    tax = ps_rela %>% vegan_tax() %>%
      as.data.frame()
    head(tax)
    tax$ID = row.names(tax)
    p.row.lab= ggplot(tax, aes(y = ID, x=1, fill=!!sym(row.lab))) + geom_tile() +
      scale_fill_brewer(palette = 'Set3',name="Cell Type") +
      theme_void()
  }



  map = phyloseq::sample_data(ps_rela) %>% as.tibble() %>%
    dplyr::arrange(Group) %>% as.data.frame()
  map$ID
  pcm = reshape2::melt(data, id = c("id"))
  pcm$variable = factor(pcm$variable,levels = map$ID)

  if (ord.col == TRUE) {
    pcm$variable = factor(pcm$variable,levels = axis_order.s)

  }

  pcm$id = factor(pcm$id,levels = rig$id)
  head(pcm)
  col0 =colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral")[11:1])(60)
  col1 = ggsci::pal_gsea()(12)

  p1 = ggplot(pcm, aes(y = id, x = variable)) +
    # geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) +
    geom_tile(aes(size = value,fill = value))+
    scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) +
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  +
    # scale_fill_manual(values = colours, guide = FALSE) +
    scale_x_discrete(limits = rev(levels(pcm$variable)))  +
    scale_y_discrete(position = "right") +
    scale_fill_gradientn(colours =col1)+
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
    scale_fill_gradientn(colours = col1) +
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.x = element_text(colour = "black",angle = 90)

    )

  p1 <- p1  %>%
    aplot::insert_right(p_rig, width=.2)
  p2 <- p2  %>%
    aplot::insert_right(p_rig, width=.2)

  if (!is.null(row.lab)) {
    p1 <- p1  %>%
      aplot::insert_left(p.row.lab, width=.02)
    p2 <- p2  %>%
      aplot::insert_left(p.row.lab, width=.02)
  }

  if (row_cluster ==  TRUE) {
    p1 <- p1  %>%
      aplot::insert_left(ggtree_plot, width=.2)
    p2 <- p2  %>%
      aplot::insert_left(ggtree_plot, width=.2)
  }






  if (label ==  TRUE) {
    p1 <- p1  %>%
      aplot::insert_top(labels, height=.02)
    p2 <- p2  %>%
      aplot::insert_top(labels, height=.02)
  }

  if (col_cluster ==  TRUE) {
    p1 <- p1  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
    p2 <- p2  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
  }

  return(list(p1,p2))
}




Microheatmap3 <- function(ps_rela,
                         id,
                         label =  TRUE,
                         col_cluster = TRUE,
                         row_cluster = TRUE,
                         ord.col = TRUE,
                         axis_order.s = axis_order.s,
                         row.lab = NULL
){

  map = phyloseq::sample_data(ps_rela)
  map$ID = row.names(map)
  phyloseq::sample_data(ps_rela) = map
  otu = as.data.frame(t(ggClusterNet::vegan_otu(ps_rela)))
  otu = as.matrix(otu[id,])
  ps_heatm = phyloseq::phyloseq(
    phyloseq::otu_table(otu,taxa_are_rows = TRUE),
    phyloseq::tax_table(ps_rela),
    phyloseq::sample_data(ps_rela)

  )

  # print(ps_heatm)
  datah <- as.data.frame(t(ggClusterNet::vegan_otu(ps_heatm)))
  head(datah)
  tax = as.data.frame(ggClusterNet::vegan_tax(ps_heatm))

  otutaxh = cbind(datah,tax)
  head(otutaxh)

  otutaxh$id = paste(row.names(otutaxh))
  # otutaxh$id  = row.names(otutaxh)
  row.names(otutaxh) = otutaxh$id


  data <- otutaxh[,c("id",phyloseq::sample_names(ps_rela))]

  rig <- data[,phyloseq::sample_names(ps_rela)] %>% rowMeans() %>% as.data.frame()
  head(rig)
  colnames(rig) = "MeanAbundance"
  rig$id = row.names(rig)
  rig = rig %>% dplyr::arrange(MeanAbundance)
  rig$id = factor(rig$id,levels = rig$id)
  p_rig = ggplot(rig) + geom_bar(aes(y = id,x = MeanAbundance),
                                 fill = "#A54657",
                                 stat = "identity") + theme_void()

  tem = data[,phyloseq::sample_names(ps_rela)] %>% as.matrix()

  tem = scale(t(tem)) %>% t() %>%
    as.data.frame()
  data[,phyloseq::sample_names(ps_rela)] = tem


  # data[data > 0.3]<-0.3
  mat <- data[,-1] #drop gene column as now in rows

  if (row_cluster ==  TRUE) {
    clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix
    ggtree_plot <- ggtree::ggtree(clust)
  }
  if (col_cluster ==  TRUE) {
    v_clust <- hclust(dist(mat %>% as.matrix() %>% t()))
    ggtree_plot_col <- ggtree::ggtree(v_clust) + ggtree::layout_dendrogram()
  }

  if (label ==  TRUE) {
    map = as.data.frame(phyloseq::sample_data(ps_rela))
    map$ID = row.names(map)
    labels= ggplot(map, aes(x = ID, y=1, fill=Group)) + geom_tile() +
      scale_fill_brewer(palette = 'Set1',name="Cell Type") +
      theme_void()
  }

  if (!is.null(row.lab) ) {
    tax = ps_rela %>% vegan_tax() %>%
      as.data.frame()
    head(tax)
    tax$ID = row.names(tax)
    p.row.lab= ggplot(tax, aes(y = ID, x=1, fill=!!sym(row.lab))) + geom_tile() +
      scale_fill_brewer(palette = 'Set2',name="Cell Type") +
      theme_void()
  }



  map = phyloseq::sample_data(ps_rela) %>% as.tibble() %>%
    dplyr::arrange(Group) %>% as.data.frame()
  map$ID

  tax = ps_rela %>% vegan_tax() %>% as.data.frame()
  data$id =  tax[,new.id]
  pcm = reshape2::melt(data, id = c("id"))
  pcm$variable = factor(pcm$variable,levels = map$ID)

  if (ord.col == TRUE) {
    pcm$variable = factor(pcm$variable,levels = axis_order.s)

  }

  # pcm$id = factor(pcm$id,levels = rig$id)
  head(pcm)
  col0 =colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral")[11:1])(60)
  col1 = ggsci::pal_gsea()(12)

  p1 = ggplot(pcm, aes(y = id, x = variable)) +
    # geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) +
    geom_tile(aes(size = value,fill = value))+
    scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) +
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  +
    # scale_fill_manual(values = colours, guide = FALSE) +
    scale_x_discrete(limits = rev(levels(pcm$variable)))  +
    scale_y_discrete(position = "right") +
    scale_fill_gradientn(colours =col1)+
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
    scale_fill_gradientn(colours = col1) +
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.x = element_text(colour = "black",angle = 90)

    )

  p1 <- p1  %>%
    aplot::insert_right(p_rig, width=.2)
  p2 <- p2  %>%
    aplot::insert_right(p_rig, width=.2)

  if (!is.null(row.lab)) {
    p1 <- p1  %>%
      aplot::insert_left(p.row.lab, width=.02)
    p2 <- p2  %>%
      aplot::insert_left(p.row.lab, width=.02)
  }

  if (row_cluster ==  TRUE) {
    p1 <- p1  %>%
      aplot::insert_left(ggtree_plot, width=.2)
    p2 <- p2  %>%
      aplot::insert_left(ggtree_plot, width=.2)
  }






  if (label ==  TRUE) {
    p1 <- p1  %>%
      aplot::insert_top(labels, height=.02)
    p2 <- p2  %>%
      aplot::insert_top(labels, height=.02)
  }

  if (col_cluster ==  TRUE) {
    p1 <- p1  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
    p2 <- p2  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
  }

  return(list(p1,p2))
}




Microheatmap.upper <- function(ps_rela,
                         id,
                         id.group = id.ko.g,
                         label =  TRUE,
                         label.r = TRUE,
                         col_cluster = TRUE,
                         row_cluster = TRUE,
                         ord.col = TRUE,
                         axis_order.s = axis_order.s,
                         row.lab = NULL,
                         col1 = ggsci::pal_gsea()(12)
){

  map = phyloseq::sample_data(ps_rela)
  map$ID = row.names(map)
  phyloseq::sample_data(ps_rela) = map
  otu = as.data.frame(t(ggClusterNet::vegan_otu(ps_rela)))
  otu = as.matrix(otu[id,])
  ps_heatm = phyloseq::phyloseq(
    phyloseq::otu_table(otu,taxa_are_rows = TRUE),
    phyloseq::tax_table(ps_rela),
    phyloseq::sample_data(ps_rela)

  )

  # print(ps_heatm)
  datah <- as.data.frame(t(ggClusterNet::vegan_otu(ps_heatm)))
  head(datah)
  tax = as.data.frame(ggClusterNet::vegan_tax(ps_heatm))

  otutaxh = cbind(datah,tax)
  head(otutaxh)

  otutaxh$id = paste(row.names(otutaxh))
  # otutaxh$id  = row.names(otutaxh)
  row.names(otutaxh) = otutaxh$id


  data <- otutaxh[,c("id",phyloseq::sample_names(ps_rela))]

  rig <- data[,phyloseq::sample_names(ps_rela)] %>% rowMeans() %>% as.data.frame()
  head(rig)
  colnames(rig) = "MeanAbundance"
  rig$id = row.names(rig)
  rig = rig %>% dplyr::arrange(MeanAbundance)
  rig$id = factor(rig$id,levels = rig$id)
  p_rig = ggplot(rig) + geom_bar(aes(y = id,x = MeanAbundance),
                                 fill = "#A54657",
                                 stat = "identity") + theme_void()

  tem = data[,phyloseq::sample_names(ps_rela)] %>% as.matrix()

  tem = scale(t(tem)) %>% t() %>%
    as.data.frame()
  data[,phyloseq::sample_names(ps_rela)] = tem


  # data[data > 0.3]<-0.3
  mat <- data[,-1] #drop gene column as now in rows

  if (row_cluster ==  TRUE) {
    clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix
    ggtree_plot <- ggtree::ggtree(clust)
  }
  if (col_cluster ==  TRUE) {
    v_clust <- hclust(dist(mat %>% as.matrix() %>% t()))
    ggtree_plot_col <- ggtree::ggtree(v_clust) + ggtree::layout_dendrogram()
  }

  if (label ==  TRUE) {
    map = as.data.frame(phyloseq::sample_data(ps_rela))
    map$ID = row.names(map)
    labels= ggplot(map, aes(x = ID, y=1, fill=Group)) + geom_tile() +
      scale_fill_brewer(palette = 'Set1',name="Cell Type") +
      theme_void()
  }
  if (label.r == TRUE) {
    tem = data.frame(ID = id,Group = id.group) %>% arrange(desc(Group))
    tem$ID = factor(tem$ID,levels = tem$ID)
    labels.g.r = ggplot(tem, aes(y = ID, x=1, fill=Group)) + geom_tile() +
      scale_fill_brewer(palette = 'Set3',name="Cell Type") +
      theme_void()
  }

  if (!is.null(row.lab) ) {
    tax = ps_rela %>% vegan_tax() %>%
      as.data.frame()
    head(tax)
    tax$ID = row.names(tax)
    p.row.lab= ggplot(tax, aes(y = ID, x=1, fill=!!sym(row.lab))) + geom_tile() +
      scale_fill_brewer(palette = 'Set3',name="Cell Type") +
      theme_void()
  }



  map = phyloseq::sample_data(ps_rela) %>% as.tibble() %>%
    dplyr::arrange(Group) %>% as.data.frame()
  map$ID

  # tax = ps_rela %>% vegan_tax() %>% as.data.frame()
  # data$id =  tax[,new.id]
  pcm = reshape2::melt(data, id = c("id"))
  pcm$variable = factor(pcm$variable,levels = map$ID)

  if (ord.col == TRUE) {
    pcm$variable = factor(pcm$variable,levels = axis_order.s)

  }

  pcm$id = factor(pcm$id,levels = rig$id)
  head(pcm)

  pcm$id = factor(pcm$id,levels = tem$ID)
  # col0 =colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral")[11:1])(60)
  # col1 = ggsci::pal_gsea()(12)

  p1 = ggplot(pcm, aes(y = id, x = variable)) +
    # geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) +
    geom_tile(aes(size = value,fill = value))+
    scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) +
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  +
    # scale_fill_manual(values = colours, guide = FALSE) +
    scale_x_discrete(limits = rev(levels(pcm$variable)))  +
    scale_y_discrete(position = "right") +
    scale_fill_gradientn(colours =col1)+
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.x = element_text(colour = "black",angle = 90,hjust = 1)

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
    scale_fill_gradientn(colours = col1) +
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.x = element_text(colour = "black",angle = 90,hjust = 1)

    )

  p1 <- p1  %>%
    aplot::insert_right(p_rig, width=.2)
  p2 <- p2  %>%
    aplot::insert_right(p_rig, width=.2)

  if (!is.null(row.lab)) {
    p1 <- p1  %>%
      aplot::insert_left(p.row.lab, width=.02)
    p2 <- p2  %>%
      aplot::insert_left(p.row.lab, width=.02)
  }

  if (row_cluster ==  TRUE) {
    p1 <- p1  %>%
      aplot::insert_left(ggtree_plot, width=.2)
    p2 <- p2  %>%
      aplot::insert_left(ggtree_plot, width=.2)
  }

  if (label.r == TRUE) {
    p1 <- p1  %>%
      aplot::insert_left(labels.g.r, width=.01)
    p2 <- p2  %>%
      aplot::insert_left(labels.g.r, width=.01)


  }




  if (label ==  TRUE) {
    p1 <- p1  %>%
      aplot::insert_top(labels, height=.02)
    p2 <- p2  %>%
      aplot::insert_top(labels, height=.02)
  }

  if (col_cluster ==  TRUE) {
    p1 <- p1  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
    p2 <- p2  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
  }

  return(list(p1,p2,plotdata = pcm))
}
