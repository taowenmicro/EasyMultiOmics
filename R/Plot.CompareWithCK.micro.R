
## 准备结果存储路径
# # 准备脚本
# source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\EdgerSuper.R")
# #
# j = "Phylum"
# ps = ps
# CK = "WT"
# abun = 0.1
# result = res
#
# library(parallel)
# library(tidyverse)
#
# res = Plot.CompareWithCK(ps = ps,
#     CK = "WT",
#     j = "Genus",
#     abun = 0.001,
#     cpu = 6)
#
# res[[1]]
func = function(i){
  TF  = c()
  a <- result %>%
    dplyr::select(ends_with("level")) %>%
    dplyr::filter(row_number() == i) %>%
    as.matrix() %>%
    as.vector() %>%
    unique()
  if (length(a) ==1) {
    if (a == "nosig") {
      TF = 0
    }else {TF[i] = 1}
  } else {
    TF = 1
  }
  return(TF)
}

# x <- 1: nrow(result)
# cl <- makeCluster(cpu) # 初始化六核心集群
# clusterEvalQ(cl,library(tidyverse))
# clusterExport(cl,'result')
# results <- parLapply(cl,x,func) # lapply的并行版本
# res.df <- do.call('c',results) # 整合结果
# stopCluster(cl) # 关闭集群




Plot.CompareWithCK.micro <- function(
  ps = ps,
  CK = "OE",
  j = "Genus",
  abun = 0.001,
  cpu = 6
  ){

  res = EdgerSuper.micro(ps = ps,group  = "Group",artGroup = NULL,
                      j = j)

  result= res [[2]]
  map = sample_data(ps)
  x <- 1: nrow(result)
  results <- lapply(x,func)
  res.df <- do.call('c',results) # 整合结果

  result$TF = res.df
  # 计算微生物的平均丰度
  result$mean <- result %>%
    # filter(TF == 1) %>%
    dplyr::select(one_of(unique(map$Group))) %>%
    rowMeans() %>%
    as.vector()

  #  筛选用于展示到图上的微生物

  if (length(result$mean[result$mean > abun]) == 0) {
    Sresult <- result %>%
      dplyr::filter(TF == 1) %>%
      head(50)
  } else{

    Sresult <- result %>%
      dplyr::filter(TF == 1) %>%
      dplyr::filter(mean > abun)
  }

  Sresult$ID = row.names(Sresult)
  b <- colnames(
    result %>%
      dplyr::select(ends_with("level"))
  )

  longda <- reshape2::melt(Sresult,
                           id.vars = c("ID",b),#需要保留不参与聚合的变量,
                           measure.vars = c(as.character(unique(map$Group))),#用于聚合的变量,
                           variable.name='treat',
                           value.name='abundance') %>%
    filter(treat != "WT")

print("1")
  level = c()
  for (i in 1:nrow(longda)) {
    level[i] <- longda[i,] %>%
      dplyr::select(contains(CK)) %>%
      dplyr::select(contains(c(as.character(longda$treat[i])))) %>%
      as.matrix() %>%
      as.vector()
  }

  longda$level = level

  ck <- reshape2::melt(Sresult,
                       id.vars = c("ID",b),#需要保留不参与聚合的变量,
                       measure.vars = c(as.character(unique(map$Group))),#用于聚合的变量,
                       variable.name='treat',
                       value.name='abundance') %>%
    dplyr::filter(treat == CK)  %>%
    dplyr::select("ID","abundance")

  colnames(ck)[2] = paste("CK",colnames(ck)[2],sep = "_")

  # head(plotda)

  plotda <- longda %>% dplyr::left_join(ck)  %>%
    dplyr::mutate(level2 = abundance - CK_abundance,.keep = "all") %>%
    arrange(ID)

  plotda$level2 <- plotda$level2 > 0
  plotda$abundance[plotda$level2 == F] = -plotda$abundance[plotda$level2 == F]



  Taxonomies_x = plyr::ddply(plotda,"ID", summarize,
                             label_sd = cumsum(abundance),
                             label_y = cumsum(abundance) - 0.5*abundance)

  plotdata <- cbind(plotda,Taxonomies_x[,-1])
  head(plotdata)


  plotdata$treat = factor(plotdata$treat,levels = as.character(unique(plotdata$treat)[length(plotdata$treat):1]))

  c = c()
  for (i in 1:nrow(plotdata)) {
    if (plotdata$level[i] %in% c("enriched","depleted") ) {
      c[i] = "*"
    }
    if (plotdata$level[i] == "nosig") {
      c[i] = ""
    }
  }
  plotdata$level3 = c
  plotdata$ID = factor(plotdata$ID,levels = unique(plotdata$ID)[length( unique(plotdata$ID)):1])
  head(plotdata)
  plotdata$treat
  p <- ggplot(plotdata) +　
    geom_bar(aes(y = ID,x = abundance,group = treat,fill = treat),stat = "identity",
             color  = "black",size = 0.5) +
    geom_vline(aes(xintercept=0), colour="black") +
    geom_text(aes(y = ID,x = label_y,label = level3),color = "white") +
    labs(title = paste0(CK,"_Abundance"),y = "ASV of microbiome",
         x = paste0("Relative abundance")) + theme_bw()  +
    scale_fill_manual(values = RColorBrewer::brewer.pal(9,"Set1"))

  p

  return(list(p,plotdata))
}


