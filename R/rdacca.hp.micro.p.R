#' @title Hierarchical Partitioning with Permutations for RDA, CCA, and db-RDA
#'
#' @description
#' The `rdacca.hp.micro.p` function performs hierarchical partitioning with permutation tests to evaluate the contribution of environmental variables to variance in OTU abundance data. The function supports Redundancy Analysis (RDA), Canonical Correspondence Analysis (CCA), and Distance-based Redundancy Analysis (db-RDA). It outputs bar plots and permutation results for each method.
#'
#' @param OTU A data frame or matrix of OTU abundance data, where rows represent samples and columns represent OTUs.
#' @param env A data frame of environmental variables, where rows correspond to samples and columns represent variables.
#' @param cca A logical value indicating whether to perform hierarchical partitioning using Canonical Correspondence Analysis (CCA). Default is `FALSE`.
#' @param dbRDA A logical value indicating whether to perform hierarchical partitioning using Distance-based Redundancy Analysis (db-RDA). Default is `FALSE`.
#' @param rep An integer specifying the number of permutations for the permutation test. Default is `9`.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item `p`: A ggplot2 object visualizing the hierarchical partitioning results for RDA.
#'   \item `dat`: A data frame containing the hierarchical partitioning results for RDA.
#'   \item `p1`: A ggplot2 object visualizing the hierarchical partitioning results for CCA (if `cca = TRUE`).
#'    \item `p2`: A ggplot2 object visualizing the hierarchical partitioning results for db-RDA (if `dbRDA = TRUE`).
#' }
#' @details
#' The function performs hierarchical partitioning for variance decomposition and applies permutation tests to estimate the significance of contributions from each environmental variable. The function supports three methods:
#' \enumerate{
#'   \item **RDA (Redundancy Analysis)**:
#'     - Applies Hellinger transformation to OTU data.
#'     - Performs permutation-based hierarchical partitioning using `permu.hp`.
#'     - Generates a bar plot showing the contributions of each environmental variable.
#'   \item **CCA (Canonical Correspondence Analysis)** (if `cca = TRUE`):
#'     - Uses the original OTU abundance data without transformation.
#'     - Performs hierarchical partitioning with permutation tests.
#'     - Generates a bar plot showing the contributions of each environmental variable.
#'   \item **db-RDA (Distance-based Redundancy Analysis)** (if `dbRDA = TRUE`):
#'     - Computes Bray-Curtis distance from the OTU abundance data.
#'     - Performs hierarchical partitioning with permutation tests.
#'     - Generates a bar plot showing the contributions of each environmental variable.
#' }
#' @examples
#' \dontrun{
#' res = rdacca.hp.micro(OTU = ps.tem %>% filter_OTU_ps(500) %>%vegan_otu() %>% as.data.frame(), env = ftab[,1:5], cca = FALSE)
#' }
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @export

rdacca.hp.micro.p <- function(OTU = mite,
                              env = env1,
                             # hpPath = hpPath,
                              cca  = FALSE,
                              dbRDA = FALSE,
                              rep = 9
){
  set.seed(123)
  permu_hp <- permu.hp(dv =OTU, iv = env, method = 'RDA', type = 'adjR2', permutations = rep)
  permu_hp
  #简单作图
  print("1")
  permu_hp$Variables <- rownames(permu_hp)
  permu_hp$p <- unlist(lapply(as.character(permu_hp$'Pr(>I)'), function(x) unlist(strsplit(x, ' '))[2]))
  #如需输出层次分割结果
  dat = permu_hp %>% as.data.frame() %>%dplyr:: arrange(desc(Individual))
  dat$id = row.names(dat)

  dat$id = factor(dat$id,levels = dat$id)

  print("2")
  p = ggplot(dat) +
    geom_bar(aes(x = id,y = Individual ),stat = "identity",colour="black",fill="#9ACD32" ) +
    geom_text(aes(x = id,y = Individual ,label = `Pr(>I)`), vjust = -0.3) +
    theme_classic()
p
  # filename = paste(hpPath,'Mciro.rda.hp.p.csv',sep = "")
  # write.csv(dat,filename)
  # filename = paste(hpPath,'Micro.rda.hp.p.pdf',sep = "")
  # ggsave(filename,p)
  # filename = paste(hpPath,'Micro.rda.hp.p.jpg',sep = "")
  # ggsave(filename,p)

  ###典范对应分析（CCA）

  #不同于 RDA，CCA 中的物种丰度一般不做 Hellinger 转化处理，可直接使用原始丰度矩阵
  #使用层次分割在 CCA 中分解每个环境变量的解释，详情 ?rdacca.hp
  #本示例计算校正后的 R2，在这里默认基于 1000 次置换获取对 CCA 中校正后 R2 的估计

  if (cca == TRUE) {
    set.seed(123)
    permu_hp <- permu.hp(dv =OTU, iv = env, method = 'CCA', type = 'adjR2', permutations = rep)

    #简单作图
    permu_hp$Variables <- rownames(permu_hp)
    permu_hp$p <- unlist(lapply(as.character(permu_hp$'Pr(>I)'), function(x) unlist(strsplit(x, ' '))[2]))
    #如需输出层次分割结果
    dat = permu_hp %>% as.data.frame() %>% dplyr::arrange(desc(Individual))
    dat$id = row.names(dat)
    dat$id = factor(dat$id,levels = dat$id)

    p = ggplot(dat) +
      geom_bar(aes(x = id,y = Individual ),stat = "identity",colour="black",fill="#9ACD32" ) +
      geom_text(aes(x = id,y = Individual ,label = `Pr(>I)`), vjust = -0.3) +
      theme_classic()
    #
    # filename = paste(hpPath,'Mciro.cca.hp.p.csv',sep = "")
    # write.csv(dat,filename)
    # filename = paste(hpPath,'Micro.cca.hp.p.pdf',sep = "")
    # ggsave(filename,p)
    # filename = paste(hpPath,'Micro.cca.hp.p.jpg',sep = "")
    # ggsave(filename,p)
  }


  ###基于距离的冗余分析（db-RDA），或称典范主坐标分析（CAP）

  # #db-RDA 需要计算群落相异指数，下文以 Bray-curtis 距离为例
  # #关于 β 多样性或相异指数的计算：https://mp.weixin.qq.com/s/Jwcz2zOwL7y2eu5U3zhAWQ
  # mite.bray <- vegdist(mite, method = 'bray')
  # #db-RDA 分析环境因子对甲螨物种丰度的影响
  # #细节可参考前文：https://mp.weixin.qq.com/s/KIhGjTL1Tzc-QL7Z03LT_g
  # mite.cap <- dbrda(mite.bray~., mite.env)
  # summary(mite.cap)  #db-RDA 概要
  # exp_adj <- RsquareAdj(mite.cap)$adj.r.squared * mite.cap$CCA$eig/sum(mite.cap$CCA$eig)  #获取校正后的 R2
  # cap1_exp <- paste('CAP1:', round(exp_adj[1]*100, 2), '%')
  # cap2_exp <- paste('CAP2:', round(exp_adj[2]*100, 2), '%')
  # plot(mite.cap, display = c('wa', 'cn'), type = 'n', xlab = cap1_exp, ylab = cap2_exp)  #db-RDA 的简单作图，只显示样本点和环境变量
  # text(mite.cap, display = 'cn', col = 'blue', cex = 0.8)
  # points(mite.cap, display = 'wa', pch = 19, cex = 1)

  #使用层次分割在 db-RDA 中分解每个环境变量的解释，详情 ?rdacca.hp
  #本示例计算校正后的 R2，db-RDA 默认使用 Ezekiel 公式计算调校正后的 R2
  if (dbRDA == TRUE) {
    set.seed(123)
    permu_hp <- permu.hp(dv =OTU, iv = env, method = 'dbRDA', type = 'adjR2', permutations = rep)

    #简单作图
    permu_hp$Variables <- rownames(permu_hp)
    permu_hp$p <- unlist(lapply(as.character(permu_hp$'Pr(>I)'), function(x) unlist(strsplit(x, ' '))[2]))
    #如需输出层次分割结果
    dat = permu_hp %>% as.data.frame() %>%dplyr:: arrange(desc(Individual))
    dat$id = row.names(dat)
    dat$id = factor(dat$id,levels = dat$id)

    p = ggplot(dat) +
      geom_bar(aes(x = id,y = Individual ),stat = "identity",colour="black",fill="#9ACD32" ) +
      geom_text(aes(x = id,y = Individual ,label = `Pr(>I)`), vjust = -0.3) +
      theme_classic()
    # filename = paste(hpPath,'Mciro.dbRDA.hp.p.csv',sep = "")
    # write.csv(dat,filename)
    # filename = paste(hpPath,'Micro.dbRDA.hp.p.pdf',sep = "")
    # ggsave(filename,p)
    # filename = paste(hpPath,'Micro.dbRDA.hp.p.jpg',sep = "")
    # ggsave(filename,p)
  }

  return(list(p,dat))

}


