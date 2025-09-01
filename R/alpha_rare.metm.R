#' @title Plotting rarefaction curve for each sample or group
#' @description Input otutab and metadata, and manual set metadata column names.
#' ggplot2 show line plot, and/or standard error
#' @param otutab OTU/ASV table;
#' @param metadata matrix or data frame, including sampleID and groupID;
#' @param groupID column name for groupID, such as "Group".
#' @param start sampling OTU/ASV table with the start number sequence count;
#' @param step number of intervals for sampling
#' @param method method for calculate alpha diversity,including "observed", "chao1", "diversity_shannon", "diversity_inverse_simpson", "diversity_gini_simpson", "diversity_fisher",  "diversity_coverage", "evenness_camargo", "evenness_pielou", "evenness_simpson", "evenness_evar", "evenness_bulla", "dominance_dbp",  "dominance_dmn", "dominance_absolute","dominance_relative", "dominance_simpson", "dominance_core_abundance" ,  "dominance_gini", "rarity_log_modulo_skewness", "rarity_low_abundance", "rarity_noncore_abundance",  "rarity_rare_abundance"
#' @return ggplot2 object.
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' rare <- mean(phyloseq::sample_sums(ps.16s))/10
#' result = alpha_rare.metm(ps = ps.16s, group = "Group", method = "Richness", start = 100, step = rare)
#' p2_1 <- result[[1]]
#' p2_1
#' # Plot the rarefaction curve for a single sample
#' p2_1 <- result[[1]]
#' p2_1
#' # Provide a data table for convenient output
#' raretab <- result[[2]]
#' head(raretab)
#' # Display rarefaction curves grouped by categories
#' p2_2 <- result[[3]]
#' p2_2
#' # Plot rarefaction curves with standard deviations by groups
#' p2_3 <- result[[4]]
#' p2_3
#' @export

# 稀释曲线函数
alpha_rare.metm =function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = NULL,method = "Richness",
                          group = "Group", start = 100,step = 100){
  # 所需R包
  # library(vegan)
  # library(microbiome)
  # library(tidyverse)

  # 构建所需子函数
  #----抽平函数，输入对象为Phyloseq对象和抽平数值#----
  phyRare = function(ps = ps,N = 3000){
    # library(phyloseq)
    #---- 提取OTU表函数 #----
    vegan_otu = function(physeq){
      OTU =  otu_table(physeq)
      if(taxa_are_rows(OTU)){
        OTU =  t(OTU)
      }
      return(as(OTU,"matrix"))
    }
    #---- 提取OTU表函数
    otb = as.data.frame(t(vegan_otu(ps)))
    otb1 = vegan::rrarefy(t(otb), N)
    ps = phyloseq(otu_table(as.matrix(otb1),taxa_are_rows = F),
                  sample_data(ps)
    )
    ps
    return(ps)
  }



  ps = inputMicro(otu,tax,map,tree,ps,group  = group)

  otu = ps %>% vegan_otu() %>% round(0)
  ps = phyloseq(
    otu_table(otu,taxa_are_rows = F),
    sample_data(ps)
  )

  #---- 全部指标#----
  all = c("observed" , "chao1"  , "diversity_inverse_simpson" , "diversity_gini_simpson",
          "diversity_shannon"   ,   "diversity_fisher"   ,  "diversity_coverage"     ,    "evenness_camargo",
          "evenness_pielou"    ,   "evenness_simpson"       ,    "evenness_evar" ,   "evenness_bulla",
          "dominance_dbp"      ,  "dominance_dmn"        ,      "dominance_absolute"   ,      "dominance_relative",
          "dominance_simpson"      ,    "dominance_core_abundance" ,  "dominance_gini"  ,           "rarity_log_modulo_skewness",
          "rarity_low_abundance"   ,    "rarity_noncore_abundance",  "rarity_rare_abundance")

  #--- 运行计算#----
  for (i in seq(start,max(sample_sums(ps)), by = step) ) {

    psRe = phyRare(ps = ps, N = i)



    if (method == "Richness") {
      count = as.data.frame(t(vegan_otu(psRe)))
      # head(count)
      x = t(count) ##转置，行为样本，列为OTU
      est = vegan::estimateR(x)
      index = est[1, ]
    }

    if (method %in% c("ACE")) {
      ap_phy = estimate_richness(psRe, measures =method)
      # head(ap_phy)
      index = ap_phy$ACE
    }

    if (method %in% all) {
      alp_mic = microbiome::alpha(psRe,index=method)
      # head(alp_mic)
      index = alp_mic[,1]
    }

    tab = data.frame(ID = names(sample_sums(psRe)))
    #得到多样性的列
    tab = cbind(tab,i,index)
    # head(tab)
    if (i == start) {
      result = tab
    }
    if (i != start) {
      result = rbind(result,tab)
    }
  }

  #----稀释结果为整齐的表，#----
  for (ii in 1:length(sample_sums(ps))) {
    result$i[result$i > sample_sums(ps)[ii][[1]]]
    df_filter= filter(result, ID ==names(sample_sums(ps)[ii]) &i > sample_sums(ps)[ii][[1]])
    result$index
    result$index[result$i>sample_sums(ps)[ii][[1]]]
    a = result$i>sample_sums(ps)[ii][[1]]
    a[a == FALSE] = "a"
    b = result$ID == names(sample_sums(ps)[ii])
    b[b == FALSE] = "b"
    result$index[a== b] = NA
  }
  #----直接给列，添加分组#----
  map = as.data.frame(sample_data(ps))
  result$Group = map$Group

  ## 绘制稀释曲线
  library(ggplot2)
  main_theme =theme(panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    plot.title = element_text(vjust = -8.5,hjust = 0.1),
                    axis.title.y =element_text(size = 7,face = "bold",colour = "black"),
                    axis.title.x =element_text(size = 7,face = "bold",colour = "black"),
                    axis.text = element_text(size = 7,face = "bold"),
                    axis.text.x = element_text(colour = "black",size = 7),
                    axis.text.y = element_text(colour = "black",size = 7),
                    legend.text = element_text(size = 7,face = "bold")
  )
  p = ggplot(data= result,aes(x = i,y = index,group = ID,colour = Group)) +
    geom_smooth(span = 0.7,  method = "loess",span = 0.5,
                se = FALSE, size = 1.2) +
    labs(x= "",y=method,title="") +theme_bw()+main_theme

  #---分组求均值和标准误+se#---
  data = result
  groups= dplyr::group_by(data, Group,i)
  data2 = dplyr::summarise(groups , mean(index), sd(index))
  # head(data2)
  colnames(data2) = c(colnames(data2)[1:2],"mean","sd")
  # 按组均值绘图
  p2 = ggplot(data= data2,aes(x = i,y = mean,colour = Group)) +
    geom_smooth(span = 0.7,se = FALSE, method = "loess") +
    labs(x= "",y=method,title="") +theme_bw()+main_theme
  # p2

  # 组均值+标准差
  p4 = ggplot(data=data2,aes(x = i,y = mean,colour = Group)) +
    geom_errorbar(data = data2,aes(ymin=mean - sd, ymax=mean + sd,colour = Group),alpha = 0.4, width=.1)+labs(x= "",y=method,title="") +theme_bw()+main_theme
  p4
  return(list(p,table = result,p2,p4))
}
