
# #---三元图
# library(ggtern)
# library(tidyverse)
# library(ggClusterNet)
# library(phyloseq)
#
#
# ternpath = paste(otupath,"/ggtern/",sep = "")
# dir.create(ternpath)
# Micro_tern(ps = ps,ternpath = ternpath )
# rank.names(ps)
Micro_tern.meta = function(ps = ps

){

  ps_rela = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) );ps_rela

  otu = ggClusterNet::vegan_otu(ps_rela) %>% as.data.frame()

  #数据分组
  iris.split <- split(otu,as.factor(as.factor(phyloseq::sample_data(ps)$Group)))
  #数据分组计算平均值
  iris.apply <- lapply(iris.split,function(x)colSums(x[]))
  # 组合结果
  iris.combine <- do.call(rbind,iris.apply)
  ven2 = t(iris.combine) %>% as.data.frame()

  head(ven2)
  A <- combn(colnames(ven2),3)
  ven2$mean = rowMeans(ven2)
  tax = ggClusterNet::vegan_tax(ps)
  otutax = cbind(ven2,tax)
  head(otutax)
  colnames(otutax)[7]= "KO.1"

  otutax$Level1[otutax$Level1 == "-"] = "Unknown"
  otutax
  # i = 1

  plot = list()
  for (i in 1:dim(A)[2]) {
    x = A[1,i]
    y = A[2,i]
    z = A[3,i]
    p <- ggtern::ggtern(data=otutax,aes_string(x = x,y=y,z=z,
                                               color = "Level1",
                                               size ="mean" ))+
      geom_point() +
     theme_void()+
      theme(legend.position = "none")

    p


    plot[[i]] = p
    names(plot)[i] =  paste(paste(x,y,z,sep = "_"),sep = "")


  }
  return(list(plot,dataplot = otutax,groups = A))
}




