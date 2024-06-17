# 绘制聚类丰度图形
#
# The function named 'alpha_barplot'
# which draw barplot + error bar + test alphabet with alpha and metadata, and reture a ggplot2 object
#
# You can learn more about package at:
#
#   https://github.com/microbiota/amplicon

#' @title Plotting alpha diversity barplot for each group with anova statistics
#' @description Combine sample clustering and drawing stacked histograms of microbial categories
#' @param otu OTU/ASV table;
#' @param map Sample metadata;
#' @param tax taxonomy table
#' @param rep Number of sample replicates measured in each group
#' @param Top Select the number of microorganisms with the highest relative abundance in all samples
#' @param hcluter_method hcluter method
#' @param Group column name for groupID in map table.
#' @param cuttree cut number
#' @details
#' hcluter method is same an function hcluster
#' @return list contain ggplot object and table.
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Yong-Xin Liu \email{yxliu@@genetics.ac.cn}
#' @references
#'
#' Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai.
#' NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
#' Nature Biotechnology, 2019(37), 6:676-684, DOI: \url{https://doi.org/10.1038/s41587-019-0104-4}
#'
#' @examples
#' # data form github
#' metadata = read.table("http://210.75.224.110/github/EasyAmplicon/data/metadata.tsv", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)
#' otutab = read.table("http://210.75.224.110/github/EasyAmplicon/data/otutab.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)
#' taxonomy = read.table("http://210.75.224.110/github/EasyAmplicon/data/taxonomy.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)
#'result <-  cluMicro.bar (dist = "bray",
#'                         otu = otutab,
#'                         tax = taxonomy,
#'                         map = metadata,
#'                         rep = 6 ,# 重复数量是6个
#'                         Top = 10, # 提取丰度前十的物种注释
#'                         tran = TRUE, # 转化为相对丰度值
#'                         hcluter_method = "complete",
#'                         cuttree = 3,
#'                         Group = "Group"
#')
#'@export




cluMicro.bar.micro <- function(
        dist = "bray",
         otu = NULL,
         tax = NULL,
         map = NULL,
         tree = NULL,
         j = "Phylum", # 使用门水平绘制丰度图表
         ps = ps,
         rep = 6 ,# 重复数量是6个
         Top = 10, # 提取丰度前十的物种注释
         tran = TRUE, # 转化为相对丰度值
         hcluter_method = "complete",
         Group = "Group",
         cuttree = 3){

  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = Group)

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
  dd = merge(d,map,by = "row.names",all = F)
  row.names(dd) = dd$Row.names
  dd$Row.names = NULL


  # library(tidyverse)
  # ggtree绘图 #----
  p  = ggtree::ggtree(hc) %<+% dd +
    geom_tippoint(size=5, shape=21, aes(fill= Group, x=x)) +
    geom_tiplab(aes(color = Group,x=x * 1.2), hjust=1)
  p

  # 按照分类学门(j)合并
  psdata =  ggClusterNet::tax_glom_wt(ps = ps1_rela,ranks = j)

  # 转化丰度值
  if (tran == TRUE) {
    psdata = psdata %>% phyloseq::transform_sample_counts(function(x) {x/sum(x)} )
  }

  #--提取otu和物种注释表格
  otu = phyloseq::otu_table(psdata)
  tax = phyloseq::tax_table(psdata)


  #--按照指定的Top数量进行筛选与合并
  for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {
      tax[i,j] =tax[i,j]
    } else {
      tax[i,j]= "Other"
    }
  }
  phyloseq::tax_table(psdata)= tax

  ##转化为表格
  Taxonomies <- psdata %>%
    phyloseq::psmelt()
  head(Taxonomies)
  Taxonomies$Abundance = Taxonomies$Abundance * 100

  Taxonomies$OTU = NULL
  colnames(Taxonomies)[1] = "id"

  head(Taxonomies)

  p <- p + ggnewscale::new_scale_fill()
  p
  p1 <- facet_plot(p, panel = 'Stacked Barplot', data = Taxonomies, geom = ggstance::geom_barh,mapping = aes(x = Abundance, fill = !!sym(j)),color = "black",stat='identity' )
  p1


  grotax <- Taxonomies %>%
    dplyr::group_by(Group,!!sym(j)) %>%
    dplyr::summarise(Abundance = sum(Abundance))
  head(grotax)

  data = c()
  i = 2
  for (i in 1:length(unique(phyloseq::sample_data(psdata)$Group))) {
    a <- as.data.frame(table(phyloseq::sample_data(psdata)$Group))[i,1]

    b =  as.data.frame(table(phyloseq::sample_data(psdata)$Group))[i,2]

    c <- grotax %>%
      dplyr::filter(Group == a)
    c$Abundance <- c$Abundance/b
    head(c)
    # data = data.frame(Sample =c$Sample,Abundance = c$Abundance,aa =c$aa,Group = c$Group)
    data = c
    if (i == 1) {
      table = data
    }
    if (i != 1) {
      table = rbind(table,data)
    }

  }
 sum( grotax$Abundance)
  head(table)


#--绘制分组的聚类结果
  ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) )
  #按照分组合并OTU表格

  otu = as.data.frame((ggClusterNet::vegan_otu(ps1_rela)))




  iris.split <- split(otu,as.factor(phyloseq::sample_data(ps1_rela)$Group))
  iris.apply <- lapply(iris.split,function(x)colMeans(x))
  # 组合结果
  iris.combine <- do.call(rbind,iris.apply)
  otuG = t(iris.combine)

  ps = phyloseq::phyloseq(phyloseq::otu_table(otuG,taxa_are_rows = T),
                          phyloseq::tax_table(ps1_rela)


  )





  hc = ps %>%
    phyloseq::distance(method = dist) %>%
    stats::hclust( method = hcluter_method )


  #  take grouping with hcluster tree
  clus <- cutree(hc, cuttree )
  # 提取树中分组的标签和分组编号
  d = data.frame(label = names(clus),
                 member = factor(clus))
  # eatract mapping file

  map = data.frame(ID = unique(phyloseq::sample_data(ps1_rela)$Group),row.names = unique(phyloseq::sample_data(ps1_rela)$Group),Group = unique(phyloseq::sample_data(ps1_rela)$Group))
  # 合并树信息到样本元数据
  dd = merge(d,map,by = "row.names",all = F)
  row.names(dd) = dd$Row.names
  dd$Row.names = NULL



  # ggtree绘图 #----
  p3  = ggtree(hc) %<+% dd +
    geom_tippoint(size=5, shape=21, aes(fill= member, x=x)) +
    geom_tiplab(aes(color = member,x=x * 1.2), hjust=1)
  p3
  p3 <- p3 + ggnewscale::new_scale_fill()
  head(grotax)

  p4 <- facet_plot(p3, panel = 'Stacked Barplot', data = table, geom = ggstance::geom_barh,mapping = aes(x = Abundance, fill = !!sym(j)),color = "black",stat='identity' )
  p4

  return(list(p,p1,p3,p4,Taxonomies))
}



