#' @title Hierarchical Partitioning for RDA, CCA, and db-RDA
#'
#' @description
#' The `rdacca.hp.micro` function performs hierarchical partitioning to decompose the variance explained by environmental variables in redundancy analysis (RDA), canonical correspondence analysis (CCA), and distance-based redundancy analysis (db-RDA). It provides visualizations of individual contributions and outputs detailed hierarchical partitioning results.
#'
#' @param OTU A data frame or matrix of OTU abundance data, where rows represent samples and columns represent OTUs.
#' @param env A data frame of environmental variables, where rows correspond to samples and columns are variables.
#' @param cca A logical value indicating whether to perform hierarchical partitioning in canonical correspondence analysis (CCA). Default is `FALSE`.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item `p`: A ggplot2 object visualizing the hierarchical partitioning results for RDA.
#'   \item `dat`: A data frame containing hierarchical partitioning results for RDA or CCA.
#'   \item `p1`: A ggplot2 object visualizing the hierarchical partitioning results for db-RDA.
#'   \item `dat2`: A data frame containing hierarchical partitioning results for db-RDA.
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item **RDA (Redundancy Analysis)**:
#'     - Standardizes the OTU data using Hellinger transformation.
#'     - Conducts RDA and calculates the adjusted R² for each environmental variable using hierarchical partitioning.
#'     - Visualizes the individual contributions of environmental variables in a bar plot.
#'   \item **CCA (Canonical Correspondence Analysis)** (if `cca = TRUE`):
#'     - Uses the original OTU abundance matrix without transformation.
#'     - Performs CCA and calculates the adjusted R² for each environmental variable.
#'     - Visualizes the contributions of environmental variables in a bar plot.
#'   \item **db-RDA (Distance-based RDA)**:
#'     - Computes Bray-Curtis distances from the OTU data.
#'     - Performs db-RDA and calculates the adjusted R² for each environmental variable using hierarchical partitioning.
#'     - Visualizes the contributions of environmental variables in a bar plot.
#' }
#'
#' @examples
#' \dontrun{
#' res = rdacca.hp.micro(OTU = ps.tem %>% filter_OTU_ps(500) %>%vegan_otu() %>% as.data.frame(), env = ftab[,1:5], cca = FALSE)
#' p3_rda = res[[1]]
#' dat1 = res[[2]]
#' p3_db_rda  = res[[3]]
#' dat2 = res[[4]]
#' }
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @export
rdacca.hp.micro <- function(
                            OTU = mite,
                            env = env1,
                            #hpPath = hpPath,
                            cca  = FALSE
){
  ###冗余分析（RDA）
  mite.hel <- decostand(OTU, method = 'hellinger')
  #使用层次分割在 RDA 中分解每个环境变量的解释，详情 ?rdacca.hp
  #本示例计算校正后的 R2，RDA 默认使用 Ezekiel 公式计算调整后的 R2
  mite.rda.hp <- rdacca.hp(mite.hel, env, method = 'RDA', type = 'adjR2', scale = FALSE)
  # mite.rda.hp
  # plot(mite.rda.hp)
  #如需输出层次分割结果
  dat = mite.rda.hp$Hier.part %>% as.data.frame() %>%dplyr:: arrange(desc(Individual))
  dat$id = row.names(dat)
  dat$id = factor(dat$id,levels = dat$id)
  p = ggplot(dat) + geom_bar(aes(x = id,y = Individual ),stat = "identity",colour="black",fill="#9ACD32" ) + theme_classic()
p
  ###典范对应分析（CCA）

  #不同于 RDA，CCA 中的物种丰度一般不做 Hellinger 转化处理，可直接使用原始丰度矩阵
  #使用层次分割在 CCA 中分解每个环境变量的解释，详情 ?rdacca.hp
  #本示例计算校正后的 R2，在这里默认基于 1000 次置换获取对 CCA 中校正后 R2 的估计

  if (cca == TRUE) {
    set.seed(123)
    mite.cca.hp <- rdacca.hp(mite, env, method = 'CCA', type = 'adjR2', scale = FALSE, n.perm = 10)
    mite.cca.hp
    plot(mite.cca.hp)
    #如需输出层次分割结果
    dat = mite.cca.hp$Hier.part %>% as.data.frame() %>% arrange(desc(Individual))
    dat$id = row.names(dat)
    dat$id = factor(dat$id,levels = dat$id)
    p = ggplot(dat) + geom_bar(aes(x = id,y = Individual ),stat = "identity",colour="black",fill="#9ACD32" ) + theme_classic()
    return(list(p,dat))
  }


  ###基于距离的冗余分析（db-RDA），或称典范主坐标分析（CAP）

  # #db-RDA 需要计算群落相异指数，下文以 Bray-curtis 距离为例
  # #关于 β 多样性或相异指数的计算：https://mp.weixin.qq.com/s/Jwcz2zOwL7y2eu5U3zhAWQ
  mite.bray <- vegdist(OTU, method = 'bray')
  # #db-RDA 分析环境因子对甲螨物种丰度的影响
  # #细节可参考前文：https://mp.weixin.qq.com/s/KIhGjTL1Tzc-QL7Z03LT_g
  # mite.cap <- dbrda(mite.bray~., env)
  # summary(mite.cap)  #db-RDA 概要
  # exp_adj <- RsquareAdj(mite.cap)$adj.r.squared * mite.cap$CCA$eig/sum(mite.cap$CCA$eig)  #获取校正后的 R2
  # cap1_exp <- paste('CAP1:', round(exp_adj[1]*100, 2), '%')
  # cap2_exp <- paste('CAP2:', round(exp_adj[2]*100, 2), '%')
  # plot(mite.cap, display = c('wa', 'cn'), type = 'n', xlab = cap1_exp, ylab = cap2_exp)  #db-RDA 的简单作图，只显示样本点和环境变量
  # text(mite.cap, display = 'cn', col = 'blue', cex = 0.8)
  # points(mite.cap, display = 'wa', pch = 19, cex = 1)

  #使用层次分割在 db-RDA 中分解每个环境变量的解释，详情 ?rdacca.hp
  #本示例计算校正后的 R2，db-RDA 默认使用 Ezekiel 公式计算调校正后的 R2
  mite.cap.hp <- rdacca.hp(mite.bray, env, method = 'dbRDA', type = 'adjR2', scale = FALSE)
  #输出层次分割结果
  dat2 = mite.cap.hp$Hier.part %>% as.data.frame() %>%dplyr:: arrange(desc(Individual))
  dat2$id = row.names(dat)
  dat2$id = factor(dat$id,levels = dat$id)
  p1 = ggplot(dat) + geom_bar(aes(x = id,y = dat$`I.perc(%)`  ),stat = "identity",colour="black",fill="#9ACD32" ) + theme_classic()
  p1
  return(list(p,dat,p1,dat2))

}


