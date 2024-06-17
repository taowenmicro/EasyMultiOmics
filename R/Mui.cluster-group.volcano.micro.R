

Mui.cluster.volcano.micro = function(res = res,rs.k = 10){
  colnames(res)[match(colnames(res)[str_detect(colnames(res),"logFC")],colnames(res))] = "logFC"
  colnames(res)[match(colnames(res)[str_detect(colnames(res),"level")],colnames(res))] = "level"

  res$ID = row.names(res)
  #--对微生物进行聚类分组


  #-就算微生物矩阵
  otu = ps %>% vegan_otu() %>% t()
  rs.dist <- dist(otu)
  rs.clust <- hclust(rs.dist, method = "ward.D")
  rs.ord.names <- rs.clust$labels[rs.clust$order]
  rs.ord <- data.frame(ID = rs.ord.names)
  rs.cut <- cutree(rs.clust[c(1,2,4)], k = rs.k)
  rs.ord$Cluster <- as.factor(rs.cut[rs.ord$ID])
  head(rs.ord)

  datv = res %>% dplyr::left_join(rs.ord,by = "ID")
  head(datv)

  # for循环挑选每个cluster的top前5 gene symbol

  tm <- function(data){
    for (i in c(1:10)) {
      tem = filter(data,Cluster==i,level != "nosig") %>%
        distinct(ID,.keep_all = TRUE) %>%
        top_n(5,abs(logFC))
      if (i ==1) {
        tem2 = tem
      } else {
        tem2 = rbind(tem2,tem)
      }
    }
    return(tem2)
  }

  top <- tm(datv)
  # 先画背景柱，根据数据log2FC的max值,min值来确定
  #根据数据中log2FC区间确定背景柱长度：

  head(datv)

  tem = datv %>%
    group_by(Cluster) %>%
    dplyr::summarise(max = max(logFC),min = min(logFC)) %>%
    as.data.frame()

  col1<-data.frame(x=tem$Cluster,
                   y=tem$max)
  col2<-data.frame(x=tem$Cluster,
                   y=tem$min)
  # 绘制背景柱
  p1 <- ggplot()+
    geom_col(data = col1,
             mapping = aes(x = x,y = y),
             fill = "#dcdcdc",alpha = 0.6)+
    geom_col(data = col2,
             mapping = aes(x = x,y = y),
             fill = "#dcdcdc",alpha = 0.6)
  p1



  #把散点火山图叠加到背景柱上：
  head(datv)
  p2 <- ggplot()+
    geom_col(data = col1,
             mapping = aes(x = x,y = y),
             fill = "#dcdcdc",alpha = 0.6)+
    geom_col(data = col2,
             mapping = aes(x = x,y = y),
             fill = "#dcdcdc",alpha = 0.6)+
    geom_jitter(data = datv,
                aes(x =Cluster , y = logFC, color =level ),
                size = 1,
                width =0.4)+
    scale_color_manual(name=NULL,
                       values = c("#4393C3","#FC4E2A","grey40"))+
    labs(x="",y="log2(FoldChange)")
  p2

  # 添加X轴的分组色块标签：
  dfcol<-data.frame(x=tem$Cluster,
                    y=0,
                    label=tem$Cluster)
  # 添加分组色块标签
  dfcol$group <- c(paste("Cluster_",tem$Cluster,sep = ""))
  # 加载包
  library(RColorBrewer)
  library(MetBrewer)
  # BiocManager::install("MetBrewer")
  # 自定义分组色块的颜色
  tile_color <- met.brewer("Thomas",length(tem$Cluster))

  # 在图中镶嵌色块
  p3 <- p2 + geom_tile(data = dfcol,
                       aes(x=x,y=y),
                       height=1.75,
                       color = "black",
                       fill = tile_color,
                       alpha = 0.6,
                       show.legend = F)+
    geom_text(data=dfcol,
              aes(x=x,y=y,label=group),
              size =3.5,
              color ="white") + theme_classic()
  p3

  library(ggrepel)
  p4<-p3+geom_text_repel(
    data=top,
    aes(x=Cluster,y=logFC,label=ID),
    force = 1.2,
    arrow = arrow(length = unit(0.008, "npc"),
                  type = "open", ends = "last"))
  p4
  # 去除背景，美化图片
  p5 <- p4+
    theme_minimal()+
    theme(
      axis.title = element_text(size = 18,
                                color = "black",
                                face = "bold"),
      axis.line.y = element_line(color = "black",
                                 size = 1.2),
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.direction = "vertical",
      legend.justification = c(1,0),
      legend.text = element_text(size = 12)
    )
  p5
  return(list(p5,p3,datv,top))
}


