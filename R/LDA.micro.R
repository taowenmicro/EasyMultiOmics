



# p1 <- p_base(ps,Top = 100)
#
#
# tablda = LDA_Micro(ps = ps,Top = 100,p.lvl = 0.05,lda.lvl = 2,seed = 11, adjust.p = T)
#
# # 注释树
#
# p2 <- clade.anno_wt(p1, tablda[[1]], alpha=0.3,anno.depth = 7)
#
# ggsave("./cs6.pdf",p2,width = 15,height = 10)






LDA.micro = function(ps = ps,
                     group = "Group",
                     Top = 100,
                     p.lvl = 0.05,
                     lda.lvl = 2,
                     seed = 11,
                     adjust.p = F
                     ){

  ps = ps %>%  filter_taxa(function(x) sum(x ) > 0 , TRUE)
  alltax = ps %>%
    ggClusterNet::filter_OTU_ps(Top) %>%
    ggClusterNet::vegan_tax() %>%
    as.data.frame()
  alltax$OTU = row.names(alltax)

  alltax$Kingdom = paste(alltax$Kingdom,sep = "_Rank_")
  alltax$Phylum = paste(alltax$Kingdom,alltax$Phylum,sep = "_Rank_")
  alltax$Class = paste(alltax$Phylum,alltax$Class,sep = "_Rank_")
  alltax$Order = paste(alltax$Class,alltax$Order,sep = "_Rank_")
  alltax$Family = paste(alltax$Order,alltax$Family,sep = "_Rank_")
  alltax$Genus = paste(alltax$Family,alltax$Genus,sep = "_Rank_")
  alltax$Species = paste(alltax$Genus,alltax$Species,sep = "_Rank_")


  otu = ps %>%
    ggClusterNet::filter_OTU_ps(Top) %>%
    ggClusterNet::vegan_otu() %>%
    t() %>%
    as.data.frame()

  otu_tax = merge(otu,alltax,by = "row.names",all = F)
  head(otu_tax)

  rank1 <- otu_tax %>%
    dplyr::group_by(Kingdom) %>%
    dplyr::summarise_if(is.numeric, sum, na.rm = TRUE)
  colnames(rank1)[1] = "id"
  rank1$id = paste("k__",rank1$id,sep = "")
  rank2 <- otu_tax %>%
    dplyr::group_by(Phylum) %>%
    dplyr::summarise_if(is.numeric, sum, na.rm = TRUE)
  colnames(rank2)[1] = "id"
  rank2$id = paste("p__",rank2$id,sep = "")
  rank3 <- otu_tax %>%
    dplyr::group_by(Class) %>%
    dplyr::summarise_if(is.numeric, sum, na.rm = TRUE)
  colnames(rank3)[1] = "id"
  rank3$id = paste("c__",rank3$id,sep = "")

  rank4 <- otu_tax %>%
    dplyr::group_by(Order) %>%
    dplyr::summarise_if(is.numeric, sum, na.rm = TRUE)
  colnames(rank4)[1] = "id"
  rank4$id = paste("o__",rank4$id,sep = "")

  rank5 <- otu_tax %>%
    dplyr::group_by(Family) %>%
    dplyr::summarise_if(is.numeric, sum, na.rm = TRUE)
  colnames(rank5)[1] = "id"
  rank5$id = paste("f__",rank5$id,sep = "")

  rank6 <- otu_tax %>%
    dplyr::group_by(Genus) %>%
    dplyr::summarise_if(is.numeric, sum, na.rm = TRUE)
  colnames(rank6)[1] = "id"
  rank6$id = paste("g__",rank6$id,sep = "")

  rank7 <- otu_tax %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise_if(is.numeric, sum, na.rm = TRUE)
  colnames(rank7)[1] = "id"
  rank7$id = paste("s__",rank7$id,sep = "")

  rank8 <- otu_tax %>%
    dplyr::group_by(OTU) %>%
    dplyr::summarise_if(is.numeric, sum, na.rm = TRUE)
  colnames(rank8)[1] = "id"
  rank8$id = paste("st__",rank8$id,sep = "")

  # 合并8个分类级
  all = rbind(rank1,rank2,rank3,rank4,rank5,rank6,rank7,rank8)
  head(all)


  #--LDA排序#--------
  data1 = as.data.frame(all)
  row.names(data1) = data1$id
  data1$id = NULL

  #-构建phylose对象

  ps_G_graphlan = phyloseq::phyloseq(phyloseq::otu_table(as.matrix(data1),taxa_are_rows = TRUE),
                                     phyloseq::sample_data(ps)) %>%  filter_taxa(function(x) sum(x ) > 0, TRUE)
  ps_G_graphlan

  #----提取OTU表格

  otu = as.data.frame((ggClusterNet::vegan_otu(ps_G_graphlan)))
  otu[otu==0] <- 1
  otu = otu[ colMeans(otu) != 1]


  map = as.data.frame(phyloseq::sample_data(ps_G_graphlan))
  # otu = (otu_table)
  claslbl= map[,group] %>% as.vector() %>% .[[1]] %>% as.factor()
  # claslbl= map$Group %>% as.factor()
  set.seed(seed)
  #KW rank sum test

  rawpvalues <- apply(otu, 2, function(x) kruskal.test(x, claslbl)$p.value);
  #--得到计算后得到的p值
  ord.inx <- order(rawpvalues)
  rawpvalues <- rawpvalues[ord.inx]
  clapvalues <- p.adjust(rawpvalues, method ="fdr")

  # p.adjust
  wil_datadf <- as.data.frame(otu[,ord.inx])


  ldares <- MASS::lda(claslbl ~ .,data = wil_datadf)
  # ldares
  ldamean <- as.data.frame(t(ldares$means))
  ldamean
  class_no <<- length(unique(claslbl))
  ldamean$max <- apply(ldamean[,1:class_no],1,max);
  ldamean$min <- apply(ldamean[,1:class_no],1,min);
  #---计算LDA
  ldamean$LDAscore <- signif(log10(1+abs(ldamean$max-ldamean$min)/2),digits=3);
  head(ldamean)

  a = rep("A",length(ldamean$max))
  for (i in 1:length(ldamean$max)) {
    name =colnames(ldamean[,1:class_no])
    a[i] = name[ldamean[,1:class_no][i,] %in% ldamean$max[i]]
  }
  ldamean$class = a

  tem1 = row.names(ldamean)
  tem1 %>% as.character()
  ldamean$Pvalues <- signif(rawpvalues[match(row.names(ldamean),names(rawpvalues))],digits=5)
  ldamean$FDR <- signif(clapvalues,digits=5)
  resTable <- ldamean
  rawNms <- rownames(resTable);
  rownames(resTable) <- gsub("`", '', rawNms);


  if (adjust.p) {
    de.Num <- sum(clapvalues <= p.lvl & ldamean$LDAscore>=lda.lvl)

  } else {
    de.Num <- sum(rawpvalues <= p.lvl & ldamean$LDAscore>=lda.lvl)
  }

  if(de.Num == 0){
    current.msg <<- "No significant features were identified with given criteria.";
  }else{
    current.msg <<- paste("A total of", de.Num, "significant features with given criteria.")
  }
  print(current.msg)
  # sort by p value
  ord.inx <- order(resTable$Pvalues, resTable$LDAscore)
  resTable <- resTable[ord.inx, ,drop=FALSE]
  resTable <- resTable[,c(ncol(resTable),1:(ncol(resTable)-1))]
  resTable <- resTable[,c(ncol(resTable),1:(ncol(resTable)-1))]

  # resTable %>% tail()
  ldamean$Pvalues[is.na(ldamean$Pvalues)] = 1
  if (adjust.p) {
    taxtree = resTable[clapvalues <=p.lvl & ldamean$LDAscore>=lda.lvl,]
  } else {
    # taxtree = resTable[ldamean$Pvalues <=p.lvl & ldamean$LDAscore>=lda.lvl,]
    taxtree = resTable[ldamean$Pvalues <=p.lvl,]
  }

  #-提取所需要的颜色
  colour = c('darkgreen','red',"blue","#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
  selececol = colour[1:length(levels(as.factor(taxtree$class)))]
  names(selececol) = levels(as.factor(taxtree$class))
  A = rep("a",length(row.names(taxtree)))

  for (i in 1:length(row.names(taxtree))) {
    A[i] = selececol [taxtree$class[i]]
  }

  taxtree$color = A
  # taxtree <- taxtree[row.names(taxtree) != "k__Bacteria",]
  # node_ids <- p0$data
  # anno <- rep("white", nrow(p1$data))

  lefse_lists = data.frame(node=row.names(taxtree),
                           color=A,
                           Group = taxtree$class,
                           stringsAsFactors = FALSE
  )


  return(list(lefse_lists,taxtree))
}

# # 注释树
# p1 <- clade.anno_wt(p, lefse_lists, alpha=0.3,anno.depth = 7)
#
# ggsave("./cs7.pdf",p1,width = 15,height = 10)
#
# p1 <- clade.anno_wt(p, lefse_lists, alpha=0.3,anno.depth = 8)
# p1
#
# ggsave("./cs9.pdf",p1,width = 10,height = 10)


# gtree = p1
# anno.data = tablda
# alpha = 0.3
# anno.depth = 7
# anno.x = 10
# anno.y = 40

# gtree = p1
# anno.data = tablda[[1]]
clade.anno_wt <- function(gtree, anno.data, alpha = 0.2, anno.depth = 5, anno.x = 10,
          anno.y = 40){
  short.labs <- c(letters,paste(letters,1:500,sep = ""))
  get_offset <- function(x) {
    (x * 0.2 + 0.2)^2
  }
  get_angle <- function(node) {
    data <- gtree$data
    sp <- tidytree::offspring(data, node)$node
    sp2 <- c(sp, node)
    sp.df <- data[match(sp2, data$node), ]
    mean(range(sp.df$angle))
  }


  anno.data <- dplyr::arrange(anno.data, node)
  hilight.color <- anno.data$color
  node_list <- anno.data$node
  node_ids <- (gtree$data %>% dplyr::filter(label %in% node_list) %>%
                 dplyr::arrange(label))$node
  anno <- rep("yellow", nrow(gtree$data))

  #---添加阴影#-------
  i = 1
  for (i in 1:length(node_ids)) {
    n <- node_ids[i]
    color <- hilight.color[i]
    anno[n] <- color
    mapping <- gtree$data %>% dplyr::filter(node == n)
    nodeClass <- as.numeric(mapping$nodeDepth)
    offset <- get_offset(nodeClass)
    gtree <- gtree + geom_hilight(node = n, fill = color,
                                  alpha = alpha, extend = offset)
  }
  # gtree$layers <- rev(gtree$layers)
  # gtree <- gtree + geom_point2(aes(size = I(nodeSize)), fill = anno,
  #                              shape = 21)
  short.labs.anno <- NULL
  i = 1
  # gtree$data


  #--添加标签#--------
  for (i in 1:length(node_ids)) {
    n <- node_ids[i]
    mapping <- gtree$data %>% dplyr::filter(node == n)
    nodeClass <- as.numeric(mapping$nodeDepth)
    if (nodeClass <= anno.depth) {
      lab <- short.labs[1]
      short.labs <- short.labs[-1]
      if (is.null(short.labs.anno)) {
        short.labs.anno = data.frame(lab = lab, annot = mapping$lab2,
                                     stringsAsFactors = F)
      }else {
        short.labs.anno = rbind(short.labs.anno, c(lab,mapping$lab2))
      }
    } else {
      lab <- mapping$lab2
    }

    offset <- get_offset(nodeClass) - 0.4
    angle <- get_angle(n) + 90
    gtree <- gtree + geom_cladelabel(node = n, label = lab,

                                     angle = angle, fontsize = 1 + sqrt(nodeClass),
                                     offset = offset, barsize = NA, hjust = 0.5)
  }

  if (!is.null(short.labs.anno)) {
    anno_shapes = sapply(short.labs.anno$lab, utf8ToInt)
    stable.p <- ggpubr::ggtexttable(short.labs.anno, rows = NULL,
                                    theme = ggpubr::ttheme(
                                      colnames.style = ggpubr::colnames_style(fill = "white"),
                                      tbody.style = ggpubr::tbody_style(fill = ggpubr::get_palette("RdBu", 6))
                                    ))

  } else{
    stable.p = NULL
  }


  y = (1:length(unique(anno.data$Group)))

  pleg <- ggplot() + geom_point2(aes(
    y = y,
    x = rep(1,length(unique(anno.data$color))),fill = as.factor(1:length(unique(anno.data$Group)))
  ),pch = 21,size = 2) +
    geom_text(aes(  y = y,
                    x = rep(1,length(unique(anno.data$color))),label = unique(anno.data$Group) ),
              hjust = -1
    ) + scale_fill_manual(values = unique(anno.data$color),guide = F) +
    theme_void()
    layout <- "
  AAAAAABB#
  AAAAAABB#
  AAAAAABBC
  AAAAAABBC
  AAAAAABB#
  "
  gtree <- gtree + stable.p+ pleg + plot_layout(design = layout)


}


lefse_bar = function(taxtree = tablda[[2]]){
  # taxtree = tablda[[2]]
  taxtree$ID = row.names(taxtree)
  head(taxtree)
  taxtree$ID = gsub("_Rank_",";",taxtree$ID)
  taxtree <- taxtree %>%
    arrange(class,LDAscore)
  taxtree$ID = factor(taxtree$ID,levels=taxtree$ID)
  taxtree$class = factor(taxtree$class,levels = unique(taxtree$class))

  pbar <- ggplot(taxtree) + geom_bar(aes(y =ID, x = LDAscore,fill = class),stat = "identity") +
    scale_fill_manual(values = unique(taxtree$color))  +
    scale_x_continuous(limits = c(0,max(taxtree$LDAscore)*1.2))
  return(pbar)
}






lefse_bar2 = function(ps = pst,taxtree = tablda[[2]]){
  taxtree = tablda[[2]]
  taxtree$ID = row.names(taxtree)
  head(taxtree)
  taxtree$ID = gsub("_Rank_",";",taxtree$ID)

  taxtree$rankid = sapply(strsplit(taxtree$ID, "__"), `[`, 1)
  tem = taxtree %>% filter(rankid %in% c("st"))
  tem$ASV = sapply(strsplit(tem$ID, "st__"), `[`, 2)
  row.names(tem) = tem$Genus
  tem$ID = tem$ASV
  tem <- tem %>%
    arrange(class,LDAscore)
  tem$ID = factor(tem$ID,levels=tem$ID)
  tem$class = factor(tem$class,levels = unique(tem$class))
  head(tem)

  tax.t = ps %>% subset_taxa.wt("OTU",as.character(tem$ID)) %>% vegan_tax() %>% as.data.frame() %>%
    rownames_to_column("ID")
  head(tax.t)
  tem = tem %>% left_join(tax.t,by = c("ASV"= "ID"))
  tem$ID2 = paste0(tem$Genus,"(",tem$ID,")")
  if (unique(tem$class) %>% length() == 2) {
    # a = unique(tem$class)
    # tem2 = tem %>% filter(class == a[2]) %>% mutate(LDAscore2 = LDAscore * -1)
    # tem$LDAscore[tem$class == a[2]] = tem2$LDAscore2
    # pbar <- ggplot(tem) + geom_bar(aes(y =ID, x = LDAscore,fill = class),stat = "identity") +
    #   scale_fill_manual(values = unique(taxtree$color)) + theme_classic()
    tem$just = 1
    a = unique(tem$class) %>% as.character()
    tem2 = tem %>% filter(class == a[2]) %>% mutate(LDAscore2 = LDAscore * -1,
                                                    just2 = just *0
    )
    tem$LDAscore[tem$class == a[2]] = tem2$LDAscore2
    tem$just[tem$class == a[2]] = tem2$just2
    head(tem)

    pbar <- ggplot(tem) + geom_bar(aes(y =ID, x = LDAscore,fill = class),stat = "identity") +
      scale_fill_manual(values = unique(taxtree$color))  +
      geom_text(aes(y =ID, x =0,label = ID2),hjust = tem$just) +
      theme_classic()+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line.y = element_blank()
      )

  } else{
    pbar <- ggplot(tem) + geom_bar(aes(y =ID2, x = LDAscore,fill = class),stat = "identity") +
      scale_fill_manual(values = unique(taxtree$color)) +
      scale_x_continuous(limits = c(0,max(taxtree$LDAscore)*1.2)) +
      theme_classic()+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line.y = element_blank()
      )
  }
  return(pbar)
}


