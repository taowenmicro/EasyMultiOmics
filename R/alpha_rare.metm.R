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
  # i = 19141935
  for (i in seq(start,max(sample_sums(ps)), by = step) ) {

    # psRe = phyRare(ps = ps, N = i,cores = 6)
    psRe <- phyRare_auto(ps = ps, N = i, drop_insufficient = TRUE, rngseed = 123)



    if (method == "Richness") {
      count = as.data.frame(t(vegan_otu(psRe)))
      # head(count)
      x = t(count) ##转置，行为样本，列为OTU
      # est = vegan::estimateR(x)
      # index = est[1, ]
      index  <- vegan::specnumber(x)
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

  tem = data.frame(ID = row.names(map),Group = map$Group)
  result = result %>% left_join(tem)


  p = ggplot(data= result,aes(x = i,y = index,group = ID,colour = Group)) +
    geom_smooth(
                method = "gam",
                span = 1,
                se = FALSE, size = 1.2) +
    labs(x= "",y=method,title="") + theme_nature()

  #---分组求均值和标准误+se#---
  data = result
  groups= dplyr::group_by(data, Group,i)
  data2 = dplyr::summarise(groups , mean(index), sd(index))
  # head(data2)
  colnames(data2) = c(colnames(data2)[1:2],"mean","sd")
  # 按组均值绘图
  p2 = ggplot(data= data2,aes(x = i,y = mean,colour = Group)) +
    geom_smooth(span = 0.7,
                se = FALSE,
                method = "gam"
                ) +
    labs(x= "",y=method,title="") + theme_nature()
  # p2

  # 组均值+标准差
  p4 = ggplot(data=data2,aes(x = i,y = mean,colour = Group)) +
    geom_smooth(span = 0.7,
                se = FALSE,
                method = "gam"
    ) +
    geom_errorbar(data = data2,aes(ymin=mean - sd, ymax=mean + sd,colour = Group),alpha = 0.4, width=.3)+labs(x= "",y=method,title="") +
    theme_nature()
  p4
  return(list(p,table = result,p2,p4))
}
