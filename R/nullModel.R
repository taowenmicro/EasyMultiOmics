#' @title Null Model Analysis for Microbial Communities
#'
#' @description
#' This function performs null model analysis to compare the observed beta diversity
#' similarity to the expected similarity under various null model assumptions.
#' It provides effect sizes, standardized effect sizes (SES), and ratios for different
#' groups in the microbial community dataset.
#'
#' @param otu A data frame or matrix containing OTU (Operational Taxonomic Unit) abundance data.
#' Rows represent taxa, and columns represent samples.
#' @param tax A data frame containing taxonomic annotation for each OTU.
#' @param map A data frame containing sample metadata, including group information.
#' @param tree A phylogenetic tree object (optional).
#' @param ps A `phyloseq` object containing OTU, taxonomic, and sample data.
#' If provided, this supersedes `otu`, `tax`, `map`, and `tree`.
#' @param group A character string specifying the grouping variable in the metadata.
#' Default is `"Group"`.
#' @param dist.method A character string specifying the distance method to calculate beta diversity.
#' Default is `"bray"`.
#' @param gamma.method A character string specifying the gamma diversity method.
#' Options are `"total"` or `"group"`. Default is `"total"`.
#' @param transfer Transformation method for the data (not currently implemented). Default is `"none"`.
#' @param null.model A character string specifying the null model to use.
#' Options are `"ecosphere"`, `"ecosim"`, or `"frequency"`. Default is `"ecosphere"`.
#'
#' @return A list containing:
#' \itemize{
#'   \item A data frame summarizing the observed and permuted similarities,
#'   including effect sizes and standardized effect sizes for each group.
#'   \item  A data frame containing the ratio of null model similarity to observed similarity for each group.
#'   \item `Results of an ANOVA test on the ratios between groups (if applicable).
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Splits the data by groups and calculates observed beta diversity similarity.
#'   \item Applies the specified null model (e.g., `ecosphere`, `ecosim`) to generate permuted similarity values.
#'   \item Compares observed and permuted similarity using t-tests, effect sizes, and standardized effect sizes.
#'   \item Calculates the ratio of permuted to observed similarity.
#'   \item If multiple groups are present, performs ANOVA on the similarity ratios between groups.
#' }
#'
#' @examples
#' \dontrun{
#' result <- nullModel(ps = psphy,group="Group",dist.method =  "bray",gamma.method = "total",transfer = "none",null.model = "ecosphere")
#' }
#'
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#'
#' @export
nullModel <- function(otu = NULL,
         tax = NULL,
         map = NULL,
         tree = NULL ,
         ps = NULL,
         group  = "Group",
         dist.method =  "bray",
         gamma.method = "total",
         transfer = "none",
         null.model = "ecosphere"){
  ps = inputMicro(otu,tax,map,tree,ps,group  = group)

  map = as.data.frame(sample_data(ps))
  grp1 = unique(map$Group)
  grp=list()
  ### 制作分组列表
  for (i in 1:length(grp1)) {
    grp[[i]]=rownames(map)[which(map$Group==grp1[i])]
  }
  names(grp) = grp1

  report = c()
  dat4anova = c()
  grp4anova = c()
  report.ES = c()
  report.SES = c()
  # x=17
  otu = as.data.frame(t(vegan_otu(ps)))
  otu = as.matrix(otu)
  for(x in c(1:length(grp))){  #
    #print(paste("Group",x))
    dataCK1 = otu[,grp[[x]]]

    ##delete empty rows
    if(gamma.method == "group"){
      rsum1 = rowSums(dataCK1)
      tempCK1 = which(rsum1==0)
      if(length(tempCK1)!=0) {dataCK1 = dataCK1[-tempCK1,]}
    }
    # 分组，对一组计算距离
    beta.dist = vegdist(t(dataCK1),method = dist.method)
    # 转化为相似性距离
    similarity.ob = 1 - beta.dist
    #similarity.ob.sd = sd(1-beta.dist, na.rm=TRUE)
    # 统计有多少个OTU
    gamma = nrow(dataCK1)
    #统计每个样本的OTU数量
    alpha = colSums(dataCK1>0)

    # OTU求和
    if(gamma.method == "group"){
      occur = apply(dataCK1, MARGIN=1, FUN=sum)
    }else{
      occur = apply(otu, MARGIN=1, FUN=sum)  #otu[valid.row,]
    }
    #print(paste(similarity.ob, similarity.ob.sd))

    r = 100
    # 构建样本矩阵，空矩阵
    similarity.pm = matrix(0, nrow=ncol(dataCK1), ncol=ncol(dataCK1))
    similarity.pm = as.dist(similarity.pm)

    # i = 1
    for(i in 1:r){
      #print(i)
      # 构造OTU矩阵孔阵
      PRM1 = matrix(0, ncol= ncol(dataCK1), nrow = nrow(dataCK1))

      if(null.model == "ecosphere"){
        # j = 1
        for(j in 1:ncol(dataCK1)){
          # 提取该样本otu大于0的全部otu
          aa = dataCK1[dataCK1[,j]>0,j]
          PRM1[sample(1:gamma, alpha[j], replace=FALSE, prob=occur), j] = aa
        }



      }else if(null.model == "ecosim"){
        PRM1 = randomizeMatrix(dataCK1, null.model="independentswap")
      }else if(null.model == "frequency"){
        PRM1 = randomizeMatrix(dataCK1, null.model="frequency")
      }

      # 计算抽的的矩阵的距离
      dist_pm = vegdist(t(PRM1),method = dist.method)
      # 将距离转化相似度放到之前构建的空阵中
      similarity.pm = similarity.pm + (1- dist_pm)
    }



    similarity.pm = similarity.pm/r

    #plot(density(similarity.pm[i,]))
    normality = shapiro.test(similarity.pm)#正态性检测
    nor.p = normality$p.value
    ttest = t.test(similarity.pm, similarity.ob, alternative="two.sided", paired = TRUE, conf.level = 0.95)
    tt.p = ttest$p.value
    conf.int = ttest$conf.int
    pm.mean = mean(similarity.pm)
    pm.sd = sd(similarity.pm)

    ES = log(similarity.ob) - log(similarity.pm)
    effect.size = mean(ES)
    effect.size.sd = sd(ES)
    SES = (similarity.ob - similarity.pm)/pm.sd
    sd.effect.size = mean(SES)
    sd.effect.size.sd = sd(SES)
    ratio = 1 - similarity.pm / similarity.ob
    ratio.mean = mean(ratio)
    ratio.sd = sd(ratio)
    dat4anova = c(dat4anova, as.vector(ratio))
    grp4anova = c(grp4anova, rep(names(grp)[x], length(ratio)))

    conf.int.str = paste("[",paste(signif(conf.int,digits=3),collapse="~"),"]",sep="")
    report = rbind(report, c(mean(similarity.ob),sd(similarity.ob), pm.mean,  pm.sd, conf.int.str, nor.p, tt.p , effect.size, effect.size.sd, sd.effect.size, sd.effect.size.sd, ratio.mean, ratio.sd))
    report.ES = c(report.ES, effect.size)
    report.SES = c(report.SES, sd.effect.size)
  }



  rownames(report) = grp1
  colnames(report) = c("Mean of observed similarity", "Standard deviation of observed similarity",
                       "Mean of permutated similarity", "Standard deviation of permutated similarity",
                       "95% Conf int of perm similarity", "Normality test (p) on Perm similarity",
                       "T test on Ob and Perm similarity", "Effect size (ES)", "SD of ES",
                       "Standardized effect size (SES)", "SD of SES", "Mean of Ratio", "SD of Ratio")
  head(report)

  rep = t(report)
  head(rep)

  # 这个统计量代表不同群落之间是否有差异
  ##将零模型的统计检验结果保存到文件中。
  if (length(unique(grp4anova)) > 1) {
    aov.re = aov(dat4anova ~ grp4anova)
  } else {
    aov.re = NULL
  }



  #---------------将比例保存起来备用
  ratio = data.frame(ratio = dat4anova,group = grp4anova)

  return(list(rep,ratio,aov.re))
}
