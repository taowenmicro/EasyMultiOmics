# 基于特征表绘制alpha多样性稀疏曲线 Alpha rarefracation curve
# This function named 'alpha_rare_all',
# which draw curve by sample, or by group with standard error. The inputs are OTU table and metadata, and return a ggplot2 object and plot data.

#' @title Plotting rarefaction curve for each sample or group
#' @description Input otutab and metadata, and manual set metadata column names.
#' ggplot2 show line plot, and/or standard error
#' @param otutab OTU/ASV table;
#' @param metadata matrix or data frame, including sampleID and groupID;
#' @param groupID column name for groupID, such as "Group".
#' @param start sampling OTU/ASV table with the start number sequence count;
#' @param step number of intervals for sampling
#' @param method method for calculate alpha diversity,including "observed", "chao1", "diversity_shannon", "diversity_inverse_simpson", "diversity_gini_simpson", "diversity_fisher",  "diversity_coverage", "evenness_camargo", "evenness_pielou", "evenness_simpson", "evenness_evar", "evenness_bulla", "dominance_dbp",  "dominance_dmn", "dominance_absolute","dominance_relative", "dominance_simpson", "dominance_core_abundance" ,  "dominance_gini", "rarity_log_modulo_skewness", "rarity_low_abundance", "rarity_noncore_abundance",  "rarity_rare_abundance"
#' @details
#' By default, returns a list with the curve and plot data
#' \itemize{
#' \item{most used indices: "observed", "chao1", "diversity_shannon", "evenness_simpson"}
#' \item{other used indices: "diversity_inverse_simpson", "diversity_gini_simpson", "diversity_fisher",  "diversity_coverage", "evenness_camargo", "evenness_pielou", "evenness_evar", "evenness_bulla", "dominance_dbp",  "dominance_dmn", "dominance_absolute","dominance_relative", "dominance_simpson", "dominance_core_abundance" ,  "dominance_gini", "rarity_log_modulo_skewness", "rarity_low_abundance", "rarity_noncore_abundance",  "rarity_rare_abundance"}
#' @return ggplot2 object.
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Yong-Xin Liu \email{yxliu@@genetics.ac.cn}
#' @references
#'
#' Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo & Yang Bai.
#' A practical guide to amplicon and metagenomic analysis of microbiome data.
#' Protein Cell, 2020(41), 1-16, DOI: \url{https://doi.org/10.1007/s13238-020-00724-8}
#'
#' Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai.
#' NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
#' Nature Biotechnology, 2019(37), 6:676-684, DOI: \url{https://doi.org/10.1038/s41587-019-0104-4}
#'
#' @seealso alpha_boxplot alpha_rare
#' @examples
#' # rare plot with otutab and metadata
#' result = alpha_rare_all(otu = otutab, map = metadata, group = "Group", method = "chao1", start = 200, step = 200)
#' result[[1]]# output sample curve plot
#' result[[2]]# data table
#' result[[3]]# output group curve plot
#' result[[4]]# output group curve with CI plot
#' # rare plot with phyloseq Object
#' library(phyloseq)
#' ps = phyloseq(otu_table(otutab, taxa_are_rows=TRUE), sample_data(metadata))
#' result = alpha_rare_all(ps = ps, group = "Group", method = "chao1", start = 1000, step = 1000)
#' result[[4]]# output group curve with CI plot
#' @export



alpha.rare.line.meta =function(
    otu = NULL,
    tax = NULL,
    map = NULL,
    ps = NULL,
    method = "Richness",
    group = "Group",
    start = 100,
    step = 100){
  # library(microbiome)
  phyRare = function(ps = ps,N = 3000){

    otb = as.data.frame(t(ggClusterNet::vegan_otu(ps)))
    otb1 = vegan::rrarefy(t(otb), N)
    ps = phyloseq::phyloseq(phyloseq::otu_table(as.matrix(otb1),taxa_are_rows = F),
                            phyloseq::sample_data(ps)
    )

    return(ps)
  }

  #----phyloseq#----
  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)

  #---- 全部指标#----
  all = c("observed" , "chao1"  , "diversity_inverse_simpson" , "diversity_gini_simpson",
          "diversity_shannon"   ,   "diversity_fisher"   ,  "diversity_coverage"     ,    "evenness_camargo",
          "evenness_pielou"    ,   "evenness_simpson"       ,    "evenness_evar" ,   "evenness_bulla",
          "dominance_dbp"      ,  "dominance_dmn"        ,      "dominance_absolute"   ,      "dominance_relative",
          "dominance_simpson"      ,    "dominance_core_abundance" ,  "dominance_gini"  ,           "rarity_log_modulo_skewness",
          "rarity_low_abundance"   ,    "rarity_noncore_abundance",  "rarity_rare_abundance")

  #--- 运行计算#----
  for (i in seq(start,max(phyloseq::sample_sums(ps)), by = step) ) {
    psRe = phyRare(ps = ps, N = i)



    if (method == "Richness") {
      count = as.data.frame(t(ggClusterNet::vegan_otu(psRe)))
      # head(count)
      x = t(count) ##转置，行为样本，列为OTU
      est = vegan::estimateR(x)
      index = est[1, ]
    }

    if (method %in% c("ACE")) {
      ap_phy = phyloseq::estimate_richness(psRe, measures =method)
      # head(ap_phy)
      index = ap_phy$ACE
    }

    if (method %in% all) {
      alp_mic = microbiome::alpha(psRe,index=method)
      # head(alp_mic)
      index = alp_mic[,1]
    }

    tab = data.frame(ID = names(phyloseq::sample_sums(psRe)))
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

  #----稀释结果为整齐的表，为了对应map分组#----
  for (ii in 1:length(phyloseq::sample_sums(ps))) {
    result$i[result$i > phyloseq::sample_sums(ps)[ii][[1]]]
    df_filter= filter(result, ID ==names(phyloseq::sample_sums(ps)[ii]) &i > phyloseq::sample_sums(ps)[ii][[1]])
    result$index
    result$index[result$i>phyloseq::sample_sums(ps)[ii][[1]]]
    a = result$i>phyloseq::sample_sums(ps)[ii][[1]]
    a[a == FALSE] = "a"
    b = result$ID == names(phyloseq::sample_sums(ps)[ii])
    b[b == FALSE] = "b"
    result$index[a== b] = NA
  }
  #----直接给列，添加分组#----
  map = as.data.frame(phyloseq::sample_data(ps))
  result$Group = map$Group

  ## 绘制稀释曲线
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
  head(result)
  result = result %>% dplyr::filter(index != "NA")
  result$ID = as.factor(result$ID)

  p = ggplot(data= result,aes(x = i,y = index,group = ID,color = Group)) +
    geom_line() +
    labs(x= "",y=method,title="") +theme_bw()+main_theme
  p

  data2 = data %>% head(nrow(data)/3)
  p1 = ggplot(data= result,aes(x = i,y = index,group = ID,color = Group)) +
    geom_smooth(data = data2,
                aes(x = i,y = index,group = ID,color = Group),
                method='lm',formula = y~ x+I(x^2),se = F) +
    labs(x= "",y=method,title="") +theme_bw()+main_theme
  p1

  library(EasyStat)
  #---分组求均值和标准误+se#---

  data = result
  groups= dplyr::group_by(data, Group,i)
  data2 = dplyr::summarise(groups , mean(index), sd(index))
  # head(data2)
  colnames(data2) = c(colnames(data2)[1:2],"mean","sd")

  head(data2)
  p2 = ggplot(data= data2,aes(x = i,y = mean,colour = Group)) +
    geom_line()+
    geom_errorbar(aes(x = i,y = mean,ymin = mean - sd,ymax = mean + sd),width=0.1) +
    labs(x= "",y=method,title="") +theme_bw()+main_theme
  # p2


  # 组均值+标准差
  p4 = ggplot(data=data2,aes(x = i,y = mean,colour = Group)) +
    geom_smooth(data = data2[1:nrow(data2),],
                method='lm',formula = y~ x+I(x^2),se = F) +
    geom_errorbar(aes(x = i,y = mean,ymin = mean - sd,ymax = mean + sd),width=0.1) +
    labs(x= "",y=method,title="") +theme_bw()+main_theme
  p4

  return(list(p,table = result,p2,p4,p1))
}
