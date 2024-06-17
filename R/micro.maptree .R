



library(ggClusterNet)
library(phyloseq)
library(tidyverse)
library(ggraph)
library(data.tree)
library(igraph)


# ps = ps_lefse
# Top = 100
# labtab =  NULL
# seed = 11



micro.maptree <- function(
        ps = ps,
        Top = 200,
        labtab =  NULL,
        seed = 1
){
  set.seed(seed)
  ps_sub <-filter_OTU_ps(ps = ps, Top = Top)
  ps_sub
  #--------------------整理tax文件-----------------
  #注意我们本次制作的是OTU水平的maptree
  tax = as.data.frame(vegan_tax(ps_sub))
  # rank_names(ps)
  # colnames(tax) <- c("Kingdom",	"Phylum",	"Class",	"Order",	"Family",	"Genus","Species")[1:length(colnames(tax))]
  # colnames(tax) <- c("Kingdom",	"Phylum",	"Class",	"Order",	"Family",	"Genus","Species")[1:length(colnames(tax))]


  str(tax)
  #添加水平识别标示
  tax$Kingdom = paste("D_",tax$Kingdom,sep = "--")
  tax$Phylum = paste(tax$Kingdom,"P_",tax$Phylum,sep = "--")
  tax$Class = paste(tax$Phylum,"C_",tax$Class,sep = "--")
  tax$Order = paste(tax$Class,"O_",tax$Order,sep = "--")
  tax$Family = paste(tax$Order,"F_",tax$Family ,sep = "--")
  tax$Genus = paste(tax$Family,"G_",tax$Genus,sep = "--")
  tax$Species = paste(tax$Genus,"S_",tax$Species,sep = "--")
  tax$OTU = paste("ASV",row.names(tax),sep = "--")
  # # 对不同分类水平的tax未知物种进行区分
  tax_table(ps_sub) = as.matrix(tax)
  head(tax)
  row.names(tax) = tax$OTU
  #--------------------构造边文件--------------------
  #提取分组信息
  Desep_group <- as.character(colnames(tax))
  #按照form--to格式来排布数据
  edge_t = tax[c(Desep_group[1:2])]
  dim(edge_t)
  edge_t <-  dplyr::distinct(edge_t)
  head(edge_t)
  for (i in 2:7) {
    result = tax[c(Desep_group[i:(i+1)])]
    colnames(result) = colnames(edge_t)
    result <-  dplyr::distinct(result)
    edge_t = rbind(edge_t,result)
  }
  tail(edge_t)
  colnames(edge_t) = c("from","to")
  if (length(unique(tax$Kingdom)) >1) {
    edge_t$from = as.character(edge_t$from)
    edge_t$to = as.character(edge_t$to)
    buc = data.frame(from = c("king_up","king_up"),to =c("D_Archaea","D_Bacteria") )
    row.names(buc) = c("sp_K","sp_k1")
    deg = rbind(edge_t,buc)
  }
  deg = edge_t
  #---------------------------构造基本节点信息文件
  #构造tree用于提取分级结构
  tree <-FromDataFrameNetwork(deg)
  vertices_t  <-  data.frame(
    name = unique(c(as.character(deg$from), as.character(deg$to)))
  )
  # Then I can easily get the level of each node, and add it to the initial data frame:
  mylevels <- data.frame( name=tree$Get('name'), level=tree$Get("level") )
  vertices_t <- vertices_t %>%
    dplyr::left_join(., mylevels, by=c("name"="name"))
  #这里提取名称信息
  AA = rep("A",length(vertices_t$name))
  for (i in 1:length(vertices_t$name)) {
    AA [i] = strsplit(basename(as.character(vertices_t$name)), "--")[[i]][length(strsplit(basename(as.character(vertices_t$name)), "--")[[i]])]
  }
  vertices_t$shortName<- AA
  # #设置level为2的定义为标签，其他偶读都省略掉，避免杂乱无章
  vertices_t <- vertices_t %>%
    dplyr::mutate(new_label=ifelse(level==2, shortName, NA))
  row.names(vertices_t) = vertices_t$name


  tax_table = as.data.frame(vegan_tax(ps_sub))
  otu_table = as.data.frame(t(vegan_otu(ps_sub)))
  head(tax_table)
  #计算全部均值
  otu_table$mean = rowMeans(otu_table)
  #组合结果
  inde = merge(otu_table,tax_table,by = "row.names", all = TRUE)
  head(inde)

  data_plot <- inde %>%
    dplyr::group_by(OTU)  %>%
    dplyr::summarise_if(is.numeric,sum) %>%
    dplyr::right_join(vertices_t,by = c("OTU" = "name")) %>%
    as.data.frame()

  data_plot <- data_plot %>%
    dplyr::left_join(tax_table,by = "OTU") %>%
    dplyr::distinct(OTU, .keep_all = TRUE)
  #指定大小映射列将NA值做转换为0
  asa = data_plot$mean
  data_plot$mean[is.na(asa)] = 0
  row.names(data_plot) = data_plot$OTU
  #选择是否平方根标准化丰度
  # data_plot$mean[!is.na(asa)] = sqrt(data_plot$mean[!is.na(asa)])
    if (is.null(labtab)) {
      p <- maptree_abun(labtab = labtab,deg = deg,data_plot= data_plot)
    }else if (!is.null(labtab)) {
      p <- maptree_lab(labtab = labtab,deg = deg,data_plot= data_plot)
    }
  return(list(p,data_plot))
}


maptree_abun <- function(labtab = labtab,deg = deg,data_plot= data_plot){
  if (is.null(labtab)) {
    head(deg)
    head(data_plot)
    mygraph <- igraph::graph_from_data_frame(deg, vertices= data_plot )
    #-----------------------------------设置颜色映射参数-------------------------
    #
    data = ggraph::create_layout(mygraph, layout = 'circlepack',weight = mean, sort.by = NULL, direction = "out")
    head(data)
    tail(data,20)
    #设置颜色
    data$Phylum = as.character(data$Phylum)
    data$Phylum[is.na(data$Phylum)] = "others"
    data$Phylum = factor(data$Phylum,levels =unique(data$Phylum) )
    colbar <-length(unique(data$Phylum))

    data$name

    fil_colors = colorRampPalette(c( "white","#CBD588", "#599861", "orange","#DA5724", "#508578", "#CD9BCD",
                                     "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
                                     "#8569D5", "#5E738F","#D1A33D", "#8A7C64"))(colbar)
    names(fil_colors) = unique(data$Phylum)
    fil_colors[1] = "white"
    # names(fil_colors)[1] = NA
    ?geom_node_circle
   # V( mygraph)$color = fil_colors[match(V(mygraph)$Phylum ,names(fil_colors))]
    ggplot(data = data) +
      geom_node_circle(aes(
        x0 = x,
        y0 = y,
        r = r,
        fill = Phylum ) ) +
      scale_fill_manual(values= fil_colors ) +
      scale_color_manual( values=c("0" = "white", "1" = "black", "2" = "black", "3" = "black", "4"="black", "5"="black", "6"="black", "7"="black"),guide = FALSE ) +
      geom_node_text( aes(label=new_label), size=6,repel = TRUE) +
      # geom_node_label( aes(label=new_label), size=3,repel = TRUE) +
      theme_void()


      p = ggraph::ggraph(mygraph, layout = 'circlepack',weight = mean,sort.by = NULL, direction = "out") +
        geom_node_circle(aes(fill = as.factor(depth),color = as.factor(depth) ) ) +
        scale_fill_hue() +
        scale_color_manual( values=c("0" = "white", "1" = "black", "2" = "black", "3" = "black", "4"="black", "5"="black", "6"="black", "7"="black"),guide = FALSE ) +
        geom_node_text( aes(label=new_label), size=6,repel = TRUE) +
        # geom_node_label( aes(label=new_label), size=3,repel = TRUE) +
        theme_void()

    # p = ggraph::ggraph(mygraph, layout = 'circlepack',weight = mean,sort.by = NULL, direction = "out") +
    #   geom_node_circle(aes(fill = as.factor(Phylum),color = as.factor(depth) ) ) +
    #   scale_fill_manual(values= fil_colors ) +
    #   scale_color_manual( values=c("0" = "white", "1" = "black", "2" = "black", "3" = "black", "4"="black", "5"="black", "6"="black", "7"="black"),guide = FALSE ) +
    #   geom_node_text( aes(label=new_label), size=6,repel = TRUE) +
    #   # geom_node_label( aes(label=new_label), size=3,repel = TRUE) +
    #   theme_void()

  }
  return(p)
}


maptree_lab = function(labtab = labtab,deg = deg,data_plot= data_plot){
   if (!is.null(labtab)) {
     row.names(labtab) = labtab$ID
     head(labtab)
     inde = merge(labtab,tax_table,by = "row.names", all = F)
     head(inde)

     data_plot1 <- inde %>% select(OTU,lab) %>% right_join(data_plot,by = "OTU")
     head(data_plot1)


     #指定大小映射列将NA值做转换为0
     asa = data_plot1$mean
     data_plot1$mean[is.na(asa)] = 0
     #选择是否平方根标准化丰度
     # data_plot$mean[!is.na(asa)] = sqrt(data_plot$mean[!is.na(asa)])
     head(data_plot1)
     mygraph <- graph_from_data_frame(deg, vertices= data_plot1 )
     data = ggraph::create_layout(mygraph, layout = 'circlepack',weight = mean, sort.by = NULL, direction = "out")
     head(data)
     i = "lab"
     data[[i]] = as.character(data[[i]])
     data[[i]][is.na(data[[i]])] = "Other_tax"
     unique(data[[i]])
     data[[i]] = factor(data[[i]],levels =unique(data[[i]]) )
     levels(data[[i]])
     data[i] = data[[i]]
     colbar <-length(unique(data[[i]]))
     mi = RColorBrewer::brewer.pal(9,"Set1")
     fil_colors = mi[1:colbar ]
     names(fil_colors ) = unique(data$KOWT_level)

     p = ggraph(mygraph, layout = 'circlepack',weight = mean, sort.by = NULL, direction = "out") +
       geom_node_circle(aes(fill = lab,color = as.factor(depth) ) ) +
       scale_fill_manual(values= fil_colors ) +
       scale_color_manual( values=c("0" = "white", "1" = "black", "2" = "black", "3" = "black", "4"="black", "5"="black", "6"="black", "7"="black"),guide = F ) +
       geom_node_text( aes(label=new_label), size=6,repel = TRUE) +
       # geom_node_label( aes(label=new_label), size=3,repel = TRUE) +
       theme_void() +
       theme(plot.margin = unit(rep(0,4), "cm"))#

   }
   return(p)

 }



