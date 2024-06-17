
# library(RColorBrewer)#调色板调用包
# library(ggClusterNet)

#---对于样本比较多的，尝试另外一种堆叠柱状图
# 环状堆叠柱状图
# 结合样本之间的相关性环状热图




# otupath = paste(res1path,"/OTU_220130/",sep = "");otupath
# dir.create(otupath)
#
#
# barpath = paste(otupath,"/circle_Micro_strack_bar/",sep = "");print(barpath)
# dir.create(barpath)
#
#
#
# # colset3 <- c(brewer.pal(12,"Paired"),brewer.pal(9,"Pastel1"))
# # colset1 <- c(brewer.pal(12,"Set1"))
#
# library(ggtree)
# p2 = circle_starc_bar(
#   ps = ps,
#   Top = 15,
#   dist = "bray",
#   cuttree = 3,
#   hcluter_method = "complete")
#
# FileName2 <- paste(barpath,"/a2_","_bar",".jpg", sep = "")
# ggsave(FileName2, p2, width = 10, height =8 )
#
# FileName2 <- paste(barpath,"/a2_","_bar",".pdf", sep = "")
# ggsave(FileName2, p2, width = 10, height =8 )



cir_barplot.micro = function(
  ps = ps,
  Top = 15,
  dist = "bray",
  cuttree = 3,
  hcluter_method = "complete"
){
  # phyloseq(ps)对象标准化
  otu = ps %>%
    ggClusterNet::scale_micro() %>%
    ggClusterNet::vegan_otu() %>% t() %>%
    as.data.frame()
  # 导出OTU表
  # otu = as.data.frame(t(vegan_otu(ps1_rela)))

  #计算距离矩阵
  unif = phyloseq::distance(ps %>% ggClusterNet::scale_micro() , method = dist)
  # 聚类树，method默认为complete
  hc <- stats::hclust(unif, method = hcluter_method )
  #  take grouping with hcluster tree
  clus <- cutree(hc, cuttree )

  d = data.frame(label = names(clus),
                 member = factor(clus))
  # eatract mapping file
  map = as.data.frame(phyloseq::sample_data(ps))

  dd = merge(d,map,by = "row.names",all = F)
  row.names(dd) = dd$Row.names
  dd$Row.names = NULL
  # library(tidyverse)
  # ggtree绘图 #----
  p  = ggtree::ggtree(hc, layout='circular') %<+% dd +
    geom_tippoint(size=5, shape=21, aes(fill= Group, x=x)) +
    geom_tiplab(aes(color = Group,x=x * 1.2), hjust=1,offset=0.3) + xlim(-0.5,NA)
  p

  psdata =  ggClusterNet::tax_glom_wt(ps = ps %>% ggClusterNet::scale_micro(),ranks = "Phylum" )
  # 转化丰度值
  psdata = psdata%>% phyloseq::transform_sample_counts(function(x) {x/sum(x)} )

  #--提取otu和物种注释表格
  otu = phyloseq::otu_table(psdata)
  tax = phyloseq::tax_table(psdata)
  head(tax)
  #--按照指定的Top数量进行筛选与合并
  j = "Phylum"
  for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {
      tax[i,j] =tax[i,j]
    } else {
      tax[i,j]= "Other"
    }
  }
  phyloseq::tax_table(psdata)= tax

  Taxonomies <- psdata %>% phyloseq::psmelt()

  Taxonomies$Abundance = Taxonomies$Abundance * 100

  Taxonomies$OTU = NULL
  colnames(Taxonomies)[1] = "id"

  head(Taxonomies)

  dat2 = data.frame(id = Taxonomies$id,Abundance = Taxonomies$Abundance,Phylum = Taxonomies$Phylum)
  head(dat2)

  p2 <- p +
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(
      data=dat2,
      geom=geom_bar,
      mapping=aes(
        x=Abundance,
        y=id,
        fill= Phylum
      ),
      stat="identity",
      width = 0.4,
      orientation="y",
      offset=0.9,
      pwidth=2,
      axis.params=list(
        axis = "x",
        text.angle = -45,
        hjust = 0,
        vjust = 0.5,
        nbreak = 4
      )
    )  +theme_void()

  return(list(plot = p2,plotdata = dat2))

}

