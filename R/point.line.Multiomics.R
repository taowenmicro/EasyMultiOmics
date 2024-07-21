# res = heatmap.line.omics(ps01 = ps.16s,psG = ps.ms)
# grid.draw(res[[1]])
# res[[2]]
# res[[3]]
# res[[4]]




point.line.Multiomics = function(ps01 = ps.ms,
                                   ps02 = ps.trans,
                                   ps03 = ps.micro,
                                   ps04 = ps.16s,

                                   method.scale = "rela",
                                   method.cor = "spearman",

                                   lab.1 = "ms",
                                   lab.2 = "trans",
                                   lab.3 = "micro",
                                   lab.4 = "16s",
                                   r.threshold=0.85,
                                   p.threshold= 0.05,
                                   top = 100,
                                   cv = 50

){

  # No1#--------
  id <- ps01 %>%
    ggClusterNet::scale_micro(method = method.scale) %>%
    ggClusterNet::filter_OTU_ps(top) %>%
    ggClusterNet::vegan_otu() %>%
    t() %>% as.data.frame() %>%rowCV %>%
    sort(decreasing = TRUE) %>%
    head(cv) %>%
    names()
  id

  ps_rela = ps01 %>%
    scale_micro(method = method.scale)

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
  tax = as.data.frame(ggClusterNet::vegan_tax(ps_heatm))
  otutaxh = cbind(datah,tax)
  otutaxh$id = paste(row.names(otutaxh))
  row.names(otutaxh) = otutaxh$id
  head(otutaxh)
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
  # if (scale == TRUE) {
  tem = scale(t(tem)) %>% t() %>%
    as.data.frame()
  # } else if (scale == FALSE){
  #   tem = tem
  # }

  data[,phyloseq::sample_names(ps_rela)] = tem
  # data[data > 0.3]<-0.3
  mat <- data[,-1] #drop gene column as now in rows

  map = phyloseq::sample_data(ps_rela) %>% as.tibble() %>%
    dplyr::arrange(Group) %>% as.data.frame()
  map$ID

  # tax = ps_rela %>% vegan_tax() %>% as.data.frame()
  # data$id =  tax[,new.id]
  pcm = reshape2::melt(data, id = c("id"))
  pcm$variable = factor(pcm$variable,levels = map$ID)


  # if (ord.col == TRUE) {
  #   pcm$variable = factor(pcm$variable,levels = axis_order.s)
  #
  # }

  pcm$id = factor(pcm$id,levels = rig$id)
  head(pcm)
  pcm1 = pcm
  pcm1$group = lab.1
  # col0 =colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral")[11:1])(60)
  col1 = ggsci::pal_gsea()(12)

  p1 = ggplot(pcm1, aes(y = id, x = variable)) +
    geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) +
    # geom_tile(aes(fill = value))+
    # scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) +
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  +
    # scale_fill_manual(values = colours, guide = FALSE) +
    scale_x_discrete(limits = rev(levels(pcm$variable)))  +
    scale_y_discrete(position = "left") +
    scale_fill_gradientn(colours =col1)+
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.y = element_text(size = 5),
      axis.text.x = element_text(colour = "black",angle = 90)

    )
  p1


  # No2#--------
  ps_tem = ps02 %>%
    ggClusterNet::scale_micro(method = method.scale)

  id <- ps02 %>%
    ggClusterNet::filter_OTU_ps(top) %>%
    ggClusterNet::vegan_otu() %>%
    t() %>% as.data.frame() %>%rowCV %>%
    sort(decreasing = TRUE) %>%
    head(cv) %>%
    names()
  id

  ps_rela = ps02 %>%
    scale_micro()

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
  tax = as.data.frame(ggClusterNet::vegan_tax(ps_heatm))
  otutaxh = cbind(datah,tax)
  otutaxh$id = paste(row.names(otutaxh))
  row.names(otutaxh) = otutaxh$id
  head(otutaxh)
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
  # if (scale == TRUE) {
  tem = scale(t(tem)) %>% t() %>%
    as.data.frame()
  # } else if (scale == FALSE){
  #   tem = tem
  # }

  data[,phyloseq::sample_names(ps_rela)] = tem
  # data[data > 0.3]<-0.3
  mat <- data[,-1] #drop gene column as now in rows

  map = phyloseq::sample_data(ps_rela) %>% as.tibble() %>%
    dplyr::arrange(Group) %>% as.data.frame()

  pcm = reshape2::melt(data, id = c("id"))
  pcm$variable = factor(pcm$variable,levels = map$ID)

  pcm$id = factor(pcm$id,levels = rig$id)
  pcm2 = pcm
  pcm2$group = lab.2
  # col0 =colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral")[11:1])(60)
  col1 = ggsci::pal_gsea()(12)

  p2 = ggplot(pcm2, aes(y = id, x = variable)) +
    geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) +
    # geom_tile(aes(fill = value))+
    # scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) +
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  +
    # scale_fill_manual(values = colours, guide = FALSE) +
    scale_x_discrete(limits = rev(levels(pcm$variable)))  +
    scale_y_discrete(position = "left") +
    scale_fill_gradientn(colours =col1)+
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.y = element_text(size = 5),
      axis.text.x = element_text(colour = "black",angle = 90)

    )
  p2



  # No3#--------
  ps_tem = ps03 %>%
    ggClusterNet::scale_micro(method = method.scale)

  id <- ps03 %>%
    ggClusterNet::filter_OTU_ps(top) %>%
    ggClusterNet::vegan_otu() %>%
    t() %>% as.data.frame() %>%rowCV %>%
    sort(decreasing = TRUE) %>%
    head(cv) %>%
    names()
  id

  ps_rela = ps03 %>%
    scale_micro()

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
  tax = as.data.frame(ggClusterNet::vegan_tax(ps_heatm))
  otutaxh = cbind(datah,tax)
  otutaxh$id = paste(row.names(otutaxh))
  row.names(otutaxh) = otutaxh$id
  head(otutaxh)
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
  # if (scale == TRUE) {
  tem = scale(t(tem)) %>% t() %>%
    as.data.frame()
  # } else if (scale == FALSE){
  #   tem = tem
  # }

  data[,phyloseq::sample_names(ps_rela)] = tem
  # data[data > 0.3]<-0.3
  mat <- data[,-1] #drop gene column as now in rows

  map = phyloseq::sample_data(ps_rela) %>% as.tibble() %>%
    dplyr::arrange(Group) %>% as.data.frame()

  pcm = reshape2::melt(data, id = c("id"))
  pcm$variable = factor(pcm$variable,levels = map$ID)

  pcm$id = factor(pcm$id,levels = rig$id)
  pcm3 = pcm
  pcm3$group = lab.3
  # col0 =colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral")[11:1])(60)
  col1 = ggsci::pal_gsea()(12)

  p3 = ggplot(pcm3, aes(y = id, x = variable)) +
    geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) +
    # geom_tile(aes(fill = value))+
    # scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) +
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  +
    # scale_fill_manual(values = colours, guide = FALSE) +
    scale_x_discrete(limits = rev(levels(pcm$variable)))  +
    scale_y_discrete(position = "left") +
    scale_fill_gradientn(colours =col1)+
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.y = element_text(size = 5),
      axis.text.x = element_text(colour = "black",angle = 90)

    )
  p3
  # No4#--------
  ps_tem = ps04 %>%
    ggClusterNet::scale_micro(method = method.scale)

  id <- ps04 %>%
    ggClusterNet::filter_OTU_ps(top) %>%
    ggClusterNet::vegan_otu() %>%
    t() %>% as.data.frame() %>%rowCV %>%
    sort(decreasing = TRUE) %>%
    head(cv) %>%
    names()
  id

  ps_rela = ps04 %>%
    scale_micro()

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
  tax = as.data.frame(ggClusterNet::vegan_tax(ps_heatm))
  otutaxh = cbind(datah,tax)
  otutaxh$id = paste(row.names(otutaxh))
  row.names(otutaxh) = otutaxh$id
  head(otutaxh)
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
  # if (scale == TRUE) {
  tem = scale(t(tem)) %>% t() %>%
    as.data.frame()
  # } else if (scale == FALSE){
  #   tem = tem
  # }

  data[,phyloseq::sample_names(ps_rela)] = tem
  # data[data > 0.3]<-0.3
  mat <- data[,-1] #drop gene column as now in rows

  map = phyloseq::sample_data(ps_rela) %>% as.tibble() %>%
    dplyr::arrange(Group) %>% as.data.frame()

  pcm = reshape2::melt(data, id = c("id"))
  pcm$variable = factor(pcm$variable,levels = map$ID)

  pcm$id = factor(pcm$id,levels = rig$id)
  pcm4 = pcm
  pcm4$group = lab.4
  # col0 =colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral")[11:1])(60)
  col1 = ggsci::pal_gsea()(12)

  p4 = ggplot(pcm4, aes(y = id, x = variable)) +
    geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) +
    # geom_tile(aes(fill = value))+
    # scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) +
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  +
    # scale_fill_manual(values = colours, guide = FALSE) +
    scale_x_discrete(limits = rev(levels(pcm$variable)))  +
    scale_y_discrete(position = "left") +
    scale_fill_gradientn(colours =col1)+
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.y = element_text(size = 5),
      axis.text.x = element_text(colour = "black",angle = 90)

    )
  p4


  #  联合分面热图
  pcm0 = rbind(pcm1,pcm2,pcm3,pcm4)
  head(pcm0)
  pcm0$group %>% unique()
  pcm0$group = factor(pcm0$group,levels = c(lab.1,lab.2,lab.3,lab.4))

  # library(ggh4x)
  p5 = ggplot(pcm0, aes(y = id, x = variable)) +
    geom_point(aes(fill = value,size = value), alpha = 0.75, shape = 21) +
    # geom_tile(aes(fill = value))+
    # scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) +
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  +
    # scale_fill_manual(values = colours, guide = FALSE) +
    scale_x_discrete(limits = rev(levels(pcm$variable)))  +
    scale_y_discrete(position = "left") +
    scale_fill_gradientn(colours =col1)+
    facet_wrap(.~ group,scales="free_y",ncol = 2) +
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      legend.position = "top",
      axis.text.y = element_text(size = 5),
      axis.text.x = element_text(colour = "black",angle = 90),
      panel.spacing.x = unit(5, "cm"),
      panel.spacing.y = unit(3, "cm")

    ) +
    ggh4x::facetted_pos_scales(
      y = list(
        group %in% c(lab.2,lab.4) ~ scale_y_discrete(position = "right")
      )
    )
  p5

  #  跨域相关性计算
  # pst1 = ps01 %>%
  #   scale_micro(method = method.scale)

  id <- ps01 %>%
    ggClusterNet::scale_micro(method = method.scale) %>%
    ggClusterNet::filter_OTU_ps(top) %>%
    ggClusterNet::vegan_otu() %>%
    t() %>% as.data.frame() %>%rowCV %>%
    sort(decreasing = TRUE) %>%
    head(cv) %>%
    names()
  id
  pst1 = ps01 %>%
    scale_micro(method = method.scale) %>%
    subset_taxa.wt("OTU",id)
  map = sample_data(pst1)
  map$Group = "A"
  sample_data(pst1) = map


  id <- ps02 %>%
    ggClusterNet::filter_OTU_ps(top) %>%
    ggClusterNet::vegan_otu() %>%
    t() %>% as.data.frame() %>%rowCV %>%
    sort(decreasing = TRUE) %>%
    head(cv) %>%
    names()
  id
  pst2 = ps02 %>%
    ggClusterNet::scale_micro() %>%
    subset_taxa.wt("OTU",id)
  map = sample_data(pst2)
  map$Group = "A"
  sample_data(pst2) = map


  id <- ps03 %>%
    ggClusterNet::filter_OTU_ps(top) %>%
    ggClusterNet::vegan_otu() %>%
    t() %>% as.data.frame() %>%rowCV %>%
    sort(decreasing = TRUE) %>%
    head(cv) %>%
    names()
  id
  pst3 = ps03 %>%
    ggClusterNet::scale_micro() %>%
    subset_taxa.wt("OTU",id)
  map = sample_data(pst3)
  map$Group = "A"
  sample_data(pst3) = map

  id <- ps04 %>%
    ggClusterNet::filter_OTU_ps(top) %>%
    ggClusterNet::vegan_otu() %>%
    t() %>% as.data.frame() %>%rowCV %>%
    sort(decreasing = TRUE) %>%
    head(cv) %>%
    names()
  id
  pst4 = ps04 %>%
    ggClusterNet::scale_micro() %>%
    subset_taxa.wt("OTU",id)
  map = sample_data(pst4)
  map$Group = "A"
  sample_data(pst4) = map

  ps.all1 = merge.ps(ps1 = pst1 ,
                     ps2 = pst2 ,
                     N1 = 0,
                     N2 = 0,
                     scale = F,
                     onlygroup = TRUE,#不进行列合并，只用于区分不同域
                     dat1.lab = lab.1,
                     dat2.lab = lab.2)
  ps.all1

  ps.all2 = merge.ps(ps1 = ps.all1 ,
                     ps2 = pst3 ,
                     N1 = 0,
                     N2 = 0,
                     scale = F,
                     onlygroup = TRUE,#不进行列合并，只用于区分不同域
                     dat1.lab = NULL,
                     dat2.lab = lab.3)
  ps.all2

  ps.all = merge.ps(ps1 = ps.all2 ,
                    ps2 = pst4 ,
                    N1 = 0,
                    N2 = 0,
                    scale = F,
                    onlygroup = TRUE,#不进行列合并，只用于区分不同域
                    dat1.lab = NULL,
                    dat2.lab = lab.4)
  ps.all


  #  tax = ps.all %>% vegan_tax() %>%as.data.frame()
  #
  # otu = ps.all %>% vegan_otu() %>% t() %>%as.data.frame()
  #  row.names(otu)


  res = corBionetwork.st(
    ps.st= ps.all,# phyloseq对象
    g1 = "Group",# 分组1
    g2 = NULL,# 分组2
    g3 = NULL,# 分组3
    ord.g1 = NULL, # 排序顺序
    ord.g2 = NULL, # 排序顺序
    ord.g3 = NULL, # 排序顺序
    order = NULL, # 出图每行代表的变量
    fill = "filed",
    size = "igraph.degree",method = method.cor,
    clu_method = "cluster_fast_greedy",
    select_layout = TRUE,layout_net = "model_maptree2",
    r.threshold=r.threshold,
    p.threshold= p.threshold,
    maxnode = 5,
    N= 0,
    scale = TRUE,
    env = NULL,
    bio = TRUE,
    minsize = 4,maxsize = 14)

  p4 = res[[1]]
  dat = res[[2]]
  node = dat[[2]]


  edge = dat[[3]]
  edge$from = edge$OTU_2
  edge$to = edge$OTU_1
  head(edge)


  node4  = add.id.facet(pcm0,"group")
  head(node4)
  node4$group = as.factor(node4$group)
  head(node4)
  node4$id = paste0(node4$group,"_",node4$id)

  node4$group2 = node4$group %>% as.numeric()


  tem = (levels(pcm$variable))
  point.f = tem[1]
  point.r = tem[length(tem)]
  node5 = node4 %>%
    dplyr::filter(variable ==  point.f,group %in% c(lab.1,lab.3))
  head(node5)

  node6 = node4 %>%
    dplyr::filter(variable == point.r,group %in% c(lab.2,lab.4))
  head(node6)
  node7 = rbind(node5,node6)
  head(node7)
  node7$variable  = node7$variable  %>% as.character() %>% as.factor()

  edge1 = edge %>% left_join(node7,by = c("from"="id"))  %>%
    dplyr::rename("id.facet.from" = "id.facet" )

  edge2 = edge1 %>% left_join(node7,by = c("to"="id")) %>%
    dplyr::rename("id.facet.to" = "id.facet" )

  edge2 = edge2 %>%
    dplyr::filter(!is.na(id.facet.from)) %>%
    dplyr::filter(!is.na(id.facet.to))
  head(edge2)

  edge2$id.facet.to
  edge2$id.facet.from
  g = p5
  for (i in 1:length(edge2$id.facet.from)) {
    id.from = edge2$id.facet.from[i] %>% strsplit( "[_]") %>% sapply(`[`, 2) %>% as.numeric()
    id.to  = edge2$id.facet.to[i] %>% strsplit( "[_]") %>% sapply(`[`, 2)%>% as.numeric()
    id.facet1 = edge2$group2.x[i]
    id.facet2 = edge2$group2.y[i]
    g <- line.across.facets.network(g,
                                     from=id.facet1, to=id.facet2,
                                     from_point_id=id.from,
                                     to_point_id=id.to,
                                     gp=gpar(lty=1, alpha=0.5)
    )
  }

  return(list(g,p1,p2,p3,p4,p5))
}


# grid.draw(g)
