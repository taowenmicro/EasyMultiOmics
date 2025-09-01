#' @title Calculate and Visualize Environmental Contributions in RDA/CCA
#'
#' @description
#' The `RDA_CCA_explain_percent` function performs Redundancy Analysis (RDA) or Canonical Correspondence Analysis (CCA)
#' to evaluate the contributions of environmental variables to microbial community variation. The function calculates
#' the explained and unexplained variance and generates a bar plot to visualize the contribution of individual variables.
#'
#' @param ps A `phyloseq` object containing microbial community data, including the OTU table and sample metadata.
#' @param env.dat A data frame of environmental variables, where rows correspond to samples and columns represent variables.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item `out`: A data frame summarizing the percentage contribution of each environmental variable to the constrained variance.
#'   \item `exp`: A data frame summarizing the total explained and unexplained variance as proportions.
#'   \item `p`: A ggplot2 object visualizing the percentage contributions of environmental variables as a horizontal bar plot.
#' }
#'
#' @examples
#' \dontrun{
#' library(phyloseq)
#' library(vegan)
#' library(ggplot2)
#' library(tidyr)
#' result <- RDA_CCA_explain_percent(ps = ps, env.dat = envRDA)
#'
#' # Extract results
#' out <- result[[1]]  # Environmental variable contributions
#' exp <- result[[2]]  # Explained and unexplained variance
#' plot <- result[[3]] # Bar plot
#' }
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @export
RDA_CCA_explain_percent = function(ps = ps,env.dat = envRDA){
  #--变量相对丰度标准化编号：ps1_rela

  ps_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) );ps_rela
  otu = as.data.frame(t(ggClusterNet::vegan_otu(ps_rela)))
  mapping = as.data.frame( phyloseq::sample_data(ps_rela))


  # match env and fg datasets
  samp.fg = colnames(otu)
  env.st = vegan::decostand(env.dat, method="standardize", MARGIN=2,na.rm = TRUE)#
  samp.env= rownames(env.st)
  my.env = match(samp.fg, samp.env)
  env.st2 = na.omit(env.st[my.env, ])  # omit the NA rows if without fg data
  samp.env= rownames(env.st2)
  my.fg = match(samp.env, samp.fg)
  otu = otu[, my.fg]

  # for CCA calculation
  otu = t(otu)
  C.whole = vegan::rda(otu, env.st2)
  C.whole

  #----------------------------------------选择环境变量-----------------------------
  inf_factor = vegan::vif.cca(C.whole)
  inf_factor

  # delete varable with max inflation factor
  na_env = which(is.na(inf_factor))
  if(isTRUE(length(na_env) > "0") ){
    inf_factor = inf_factor[-na_env]
  }

  max_env = which(inf_factor == max(inf_factor,na.rm = T))
  env.st4 = env.st2

  while ( inf_factor[max_env] > 20){
    env.st4 = env.st4[,-max_env]
    C.reduced = vegan::cca(otu, env.st4)
    inf_factor = vegan::vif.cca(C.reduced)
    max_env = which(inf_factor == max(inf_factor,na.rm = T ))
  }


  output2 = inf_factor ;output2

  C.whole = vegan::rda(otu, env.st4)  ##rda(otu, env.st3)
  total.chi = C.whole$tot.chi
  ind.p = array(0,dim=c(1,ncol(env.st4)))
  for(j in 1:ncol(env.st4)){
    # j = 1
    ind.par = vegan::rda(otu, env.st4[,j], env.st4[,-j])
    ind.chi = ind.par$CCA$tot.chi
    ind.per = ind.chi/total.chi
    ind.p[j] = ind.per
  }
  ind.p

  rowname = colnames(env.st4);rowname
  out = matrix(data=NA,ncol=length(colnames(env.st4)),nrow=1);out
  out = ind.p
  rownames(out) = "percent"
  colnames(out) = rowname
  out



  #------提取解释比例
  total.chi = C.whole$tot.chi;total.chi
  total.constrained = C.whole$CCA$tot.chi ; total.constrained
  # 解释的比例
  explained.percent = (total.constrained) / total.chi;explained.percent
  # 未解释的比例
  unexplained.percent = (total.chi - total.constrained) / total.chi;unexplained.percent
  exp = data.frame(ID = c("explained.percent","unexplained.percent"),count = c(explained.percent,unexplained.percent))
  exp

  out=  as.data.frame(out)
  out$all = explained.percent

  out_long =  pivot_longer(out,cols = everything() , names_to = "Microbe",      # 转换后的变量列名
                           values_to = "Percentage")

  p= ggplot(out_long, aes(x = reorder(Microbe, Percentage), y = Percentage)) +
    geom_bar(stat = "identity", fill = "skyblue", color = "black") +
    coord_flip() +  # 水平显示
    theme_minimal() +
    labs(
      title = "Microbial Percentages",
      x = "Microbe",
      y = "Percentage"
    ) +
    theme(
      axis.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )
  p
  return(list(out,exp,p))
}

# result = RDA_CCA_explain_percent(ps = ps,env.dat = envRDA)
#
# out = result[[1]]
# wxp = result[[2]]
#
# filenamea = paste(RDApath,"each_env_exp_percent.csv",sep = "")
# write.csv(out,filenamea)
#
# filenamea = paste(RDApath,"all_index_explain_percent.csv",sep = "")
# write.csv(exp,filenamea)
