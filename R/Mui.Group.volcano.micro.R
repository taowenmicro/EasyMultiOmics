


Mui.Group.volcano.micro = function(res = res){
  res$ID = row.names(res)
  datv = res
  # for循环挑选每个cluster的top前5 gene symbol
  tm.g <- function(data){
    id = data$group %>% unique()

    for (i in 1:length(id)) {
      tem = filter(data,group==id[i],level != "nosig") %>%
        distinct(ID,.keep_all = TRUE) %>%
        top_n(5,abs(logFC))
      if (i == 1) {
        tem2 = tem
      } else {
        tem2 = rbind(tem2,tem)
      }
    }
    return(tem2)
  }

  top <- tm.g(datv)
  # 先画背景柱，根据数据log2FC的max值,min值来确定
  #根据数据中log2FC区间确定背景柱长度：

  head(datv)

  tem = datv %>% group_by(group) %>%
    dplyr::summarise(max = max(logFC),min = min(logFC)) %>%
    as.data.frame()

  col1<-data.frame(x=tem$group,
                   y=tem$max)
  col2<-data.frame(x=tem$group,
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
             fill = "#dcdcdc",alpha = 0.4)+
    geom_col(data = col2,
             mapping = aes(x = x,y = y),
             fill = "#dcdcdc",alpha = 0.4)+
    geom_jitter(data = datv,
                aes(x =group , y = logFC, color =level ),
                size = 1,
                width =0.4)+
    scale_color_manual(name=NULL,
                       values = c("#4393C3","#FC4E2A","grey40"))+
    labs(x="",y="log2(FoldChange)")
  p2

  # 添加X轴的分组色块标签：
  dfcol<-data.frame(x=tem$group,
                    y=0,
                    label=tem$group)
  # 添加分组色块标签
  dfcol$group <- tem$group
  # 加载包
  library(RColorBrewer)
  library(MetBrewer)
  # BiocManager::install("MetBrewer")
  # 自定义分组色块的颜色
  tile_color <- met.brewer("Thomas",length(tem$group))

  # 在图中镶嵌色块
  # library(ggrepel)
  p3 <- p2 + geom_tile(data = dfcol,
                       aes(x=x,y=y),
                       height=1.75,
                       color = "black",
                       fill = tile_color,
                       alpha = 0.6,
                       show.legend = F)+
    ggrepel::geom_text_repel(data=dfcol,
              aes(x=x,y=y,label=group),
              size =3.5,
              color ="white") + theme_classic()
  p3

  # library(ggrepel)
  p4<-p3+geom_text_repel(
    data=top,
    aes(x=group,y=logFC,label=ID),
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
