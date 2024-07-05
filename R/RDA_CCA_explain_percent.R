





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
  return(list(out,exp))
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
