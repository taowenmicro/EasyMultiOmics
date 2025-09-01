

#' @title FEAST Microbial traceability analysis
#' @description Input otutab, metadata or phyloseq object; ; output a table object.
#' @param otu OTU/ASV table;
#' @param map Sample metadata;
#' @param group group ID;
#' @param sinkG object of sink group
#' @param sourceG object of source group
#' @details
#' By default, input phyloseq object include metadata and otutab
#' @return  a table
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#'result = FEAST.micro(ps = ps.16s,group = "Group",sinkG = "OE",sourceG = c("WT","KO"))
#' FEAST(otu = otu,map = map,group = "Group",sinkG = "WT",sourceG = c("KO","OE"))
#' FEAST(ps =ps,group = "Group",sinkG = "WT",sourceG = c("KO","OE"))
#'
#' @export
#'

# #清空内存
# rm(list=ls())
# #导入otu表格
# otu = read.delim("../../data/otutab.txt",row.names = 1)
# #导入分组文件
# map = read.delim("../../data/metadata.tsv",row.names = 1)
# head(map)
# result = FEAST(otu = otu,map = map,group = "Group",sinkG = "WT",sourceG = c("KO","OE"))
# result
# #-案例二
# ps = readRDS("../../data/ps_liu.rds")
# data(ps)
# result = FEAST(ps =ps,group = "Group",sinkG = "WT",sourceG = c("KO","OE"))
# result


FEAST.micro = function(otu = otutab,map = metadata,ps = NULL,
                       group = "Group",sinkG = "WT",
                       sourceG = c("KO","OE"),
                       rep = NULL
){
  #-
  library(tidyverse)
  library("vegan")
  library("reshape2")
  # library(EasyMicrobiome)
  #

  # source(paste(path,"FEAST-master/FEAST_src/src.R",sep = "/"))
  # #
  # if (is.null(ps) ) {
  #   head(otu)
  #   otu = as.matrix(otu)
  #   str(otu)
  #   colnames(map) = gsub(group,"AA", colnames(map))
  #   map$Group = map$AA
  #   map$Group = as.factor(map$Group )
  #
  #   ps <- phyloseq(otu_table(otu, taxa_are_rows=TRUE),
  #                  sample_data(map)
  #   )
  # }
  #
  # if (!is.null(ps) ) {
  #   ps = ps
  #   map = as.data.frame(sample_data(ps))
  #   map = map[, group]
  #   colnames(map) = "Group"
  #   map$Group = as.factor(map$Group)
  #   sample_data(ps) = map
  # }

  otu = NULL
  tax = NULL
  map = NULL
  Group = group
  ps = inputMicro(otu,tax,map,tree,ps,group  = Group)
  # 提取分组文件#----
  metadata <- as.data.frame(sample_data(ps))
  head(metadata)
  #--提取otu表格#--------
  otus <-  as.data.frame(t(vegan_otu(ps)))
  otus <- t(as.matrix(otus))


  head(metadata)
  #--将分组文件置于首列
  metadata$id = row.names(metadata)
  metadata = as.tibble(metadata)
  metadata<- dplyr::arrange(metadata, Group)


  #----提取样本名称，后续添加标记#------
  envs <- metadata$id

  #--目标ID提取#-----
  mu = metadata$id[metadata$Group==sinkG]


  # 设置FEAST运行参数#----
  EM_iterations = 1000 #default value
  different_sources_flag = 1
  Proportions_est <- list()

  #-提取每个分组测定的重复数量#-------

  if (is.null(rep)) {
    rep = length(metadata$Group)/length(unique(metadata$Group))
  } else{
    rep = rep
  }

  it = 1
  for(it in 1:rep){
    # it = 6
    #提取sink和source对应样本的位置，列的位置，方便后续提取#----
    train.ix <- which(metadata$Group%in%sourceG&metadata$id %in% metadata$id[seq(1, length(metadata$Group),
                                                                                 rep)+(it-1)])
    test.ix <- which(metadata$Group==sinkG & metadata$id == mu[it])
    #--统计source样本数量#----
    num_sources <- length(train.ix)
    num_sources
    #-输入的是原始序列文件这里进行计算最小抽平数量#----------
    COVERAGE =  min(rowSums(otus[c(train.ix, test.ix),]))  #Can be adjusted by the user

    #提取source和sink对应的otu表格并抽平
    sources <- as.matrix(vegan::rrarefy(otus[train.ix,], COVERAGE))
    sinks <- as.matrix(vegan::rrarefy(t(as.matrix(otus[test.ix,])), COVERAGE))
    # sources = sources[1,1:10]
    #  tem1[1:4,1:4]
    if (length(sourceG) == 1) {
      tem1 <- sources %>% as.data.frame()
      sources = rbind(tem1,tem1) %>% as.matrix()
      dim(sources)
      row.names(sources) = paste(sourceG,1:2,sep = "")
    }

    #打印数量信息
    print(paste("Number of OTUs in the sink sample = ",length(which(sinks > 0))))
    print(paste("Seq depth in the sources and sink samples = ",COVERAGE))
    print(paste("The sink is:", envs[test.ix]))
    #FEAST 计算主函数#-------
    str(sources)
    FEAST_output<-FEAST(source=sources, sinks = sinks,
                        env = envs[train.ix], em_itr = EM_iterations, COVERAGE = COVERAGE)
    Proportions_est[[it]] <- FEAST_output$data_prop[,1]
    #整理结果#---
    names(Proportions_est[[it]]) <- c(as.character(envs[train.ix]), "unknown")

    if(length(Proportions_est[[it]]) < num_sources +1){
      tmp = Proportions_est[[it]]
      Proportions_est[[it]][num_sources] = NA
      Proportions_est[[it]][num_sources+1] = tmp[num_sources]
    }
    print("Source mixing proportions")
    print(Proportions_est[[it]])

  }

  went = as.data.frame(Proportions_est)
  colnames(went) = mu
  head(went)
  return(table = went)
}

x<-c("vegan", "dplyr", "ggrepel", "doParallel", "foreach" ,"mgcv", "reshape2", "ggplot2", "Rcpp", "RcppArmadillo")
lapply(x, require, character.only = TRUE)
library("vegan")
library("dplyr")
library("doParallel")
library("foreach")
library("mgcv")
library("reshape2")
library("ggplot2")
library("cowplot")
library("Rcpp")
library("RcppArmadillo")
cppFunction("arma::mat schur(arma::mat& a, arma::mat& b)
            {return(a % b); }", depends="RcppArmadillo")


"change_C"<-function(newcov, X){

  X=t(as.matrix(X))
  idx = 1:dim(X)[2]

  if(sum(X) > newcov){

    while(sum(X) > newcov){
      greaterone = X > 1
      samps = 20
      if(samps > length(X[greaterone]))
        samps = length(X[greaterone])
      changeidx = sample(idx[greaterone], samps, replace = F)
      X[changeidx] = X[changeidx] - 1
    }

  }

  if(sum(X) < newcov){

    while(sum(X) < newcov){
      greaterone = X > 1
      samps = 100
      if(samps > length(X[greaterone]))
        samps = length(X[greaterone])
      changeidx = sample(idx[greaterone], samps, replace = F)
      X[changeidx] = X[changeidx] + 1
    }

  }

  return(X)
}

# x = matrix(sinks, nrow = 1)

rarefy <- function(x,maxdepth){


  if(is.null(maxdepth)) return(x)

  if(!is.element(class(x)[1], c('matrix', 'data.frame','array')))
    x <- matrix(x,nrow=nrow(x))
  nr <- nrow(x)
  nc <- ncol(x)

  for(i in 1:nrow(x)){
    if(sum(x[i,]) > maxdepth){
      prev.warn <- options()$warn
      options(warn=-1)
      s <- sample(nc, size=maxdepth, prob=x[i,], replace=T)
      options(warn=prev.warn)
      x[i,] <- hist(s,breaks=seq(.5,nc+.5,1), plot=FALSE)$counts
    }
  }
  return(x)
}

"jsdmatrix" <- function(x){
  d <- matrix(0,nrow=nrow(x),ncol=nrow(x))
  for(i in 1:(nrow(x)-1)){
    for(j in (i+1):nrow(x)){
      d[i,j] <- jsd(x[i,], x[j,])
      d[j,i] <- d[i,j]
    }
  }
  return(d)
}

"jsd" <- function(p,q){
  m <- (p + q)/2
  return((kld(p,m) + kld(q,m))/2)
}


"h"<-function(x) {y <- x[x > 0]; -sum(y * log(y))};
"mult_JSD" <- function(p,q) {h(q %*% p) - q %*% apply(p, 1, h)}

"retrands"<-function(V){
  toret<-unlist(lapply(c(V), function(x) runif(1, x+1e-12, x+1e-09)))
  return(toret)
}

"getR2"<-function(x,y){
  return((cor(x,y))^2)
}

"E"<-function(alphas, sources){
  nums<-(sapply(1:length(alphas), function(n) Reduce("+", crossprod(as.numeric(alphas[n]),as.numeric(sources[[n]])))))
  denom<-(Reduce("+", nums))
  return(nums/denom)
}

"A"<-function(alph, XO, raos){
  tmp<-crossprod(alph, XO/raos)
  tmp<-rapply(list(tmp), f=function(x) ifelse(is.nan(x),0,x), how="replace" )
  tmp<-Reduce("+",unlist(tmp))
  return(tmp)
}

"M"<-function(alphas, sources, sink, observed){

  newalphs<-c()
  rel_sink <-sink/sum(sink)

  if(sum(sources[[1]]) > 1){

    sources <-lapply(sources, function(x) x/(sum(colSums(x))))
  }


  LOs<-lapply(sources, schur, b=rel_sink)
  BOs<-t(mapply(crossprod, x=sources, y=alphas))
  BOs<-split(BOs, seq(nrow(BOs)))
  BOs<-lapply(BOs, as.matrix)
  BOs<-lapply(BOs, t)
  num_list <- list()
  source_new <- list()


  for(i in 1:length(sources)){
    num <- c()
    denom <- c()
    num<-crossprod(alphas[i], (LOs[[i]]/(Reduce("+", BOs))))
    num<-rapply(list(num), f=function(x) ifelse(is.nan(x),0,x), how="replace" ) #replace na with zero
    num_list[[i]]<- num[[1]][1,] + observed[[i]][1,]

    denom <- Reduce("+",unlist(num_list[[i]]))
    source_new[[i]] <- num_list[[i]]/denom
    source_new[[i]][is.na(source_new[[i]])] = 0
  }

  sources = source_new

  newalphs<-c()
  #sink<-as.matrix(sink); #src1<-as.matrix(sources[[1]]); src2<-as.matrix(sources[[2]])
  sources<-lapply(sources, t)
  XOs<-lapply(sources,schur, b=rel_sink)
  AOs<-t(mapply(crossprod, x=sources, y=alphas))
  AOs<-split(AOs, seq(nrow(AOs)))
  AOs<-lapply(AOs, as.matrix)
  AOs<-lapply(AOs, t)
  newAs<-c()
  for(i in 1:length(sources)){
    newA<-crossprod(alphas[i], (XOs[[i]]/(Reduce("+", AOs))))
    newA<-rapply(list(newA), f=function(x) ifelse(is.nan(x),0,x), how="replace" )
    newA<-Reduce("+",unlist(newA))
    newAs<-c(newAs, newA)
  }
  tot<-sum(newAs)
  Results <- list (new_alpha = newAs/(tot), new_sources = sources)
  return(Results)
}

"do_EM"<-function(alphas, sources, observed, sink, iterations){

  curalphas<-alphas
  newalphas<-alphas
  m_guesses<-c(alphas[1])
  for(itr in 1:iterations){

    curalphas<-E(newalphas, sources)
    tmp <- M(alphas = curalphas, sources = sources, sink = sink, observed = observed)
    newalphas <- tmp$new_alpha
    sources <- tmp$new_sources

    m_guesses<-c(m_guesses, newalphas[1])
    if(abs(m_guesses[length(m_guesses)]-m_guesses[length(m_guesses)-1])<=10^-6) break

  }
  toret<-c(newalphas)
  results <- list(toret = toret, sources = sources)

  return(results)
}

"M_basic"<-function(alphas, sources, sink){
  newalphs<-c()
  XOs<-lapply(sources,schur, b=sink)
  AOs<-t(mapply(crossprod, x=sources, y=alphas))
  AOs<-split(AOs, seq(nrow(AOs)))
  AOs<-lapply(AOs, as.matrix)
  AOs<-lapply(AOs, t)
  newAs<-c()
  for(i in 1:length(sources)){
    newA<-crossprod(alphas[i], (XOs[[i]]/(Reduce("+", AOs))))
    newA<-rapply(list(newA), f=function(x) ifelse(is.nan(x),0,x), how="replace" )
    newA<-Reduce("+",unlist(newA))
    newAs<-c(newAs, newA)
  }
  tot<-sum(newAs)
  return(newAs/(tot))
}

"do_EM_basic"<-function(alphas, sources, sink, iterations){
  curalphas<-alphas
  newalphas<-alphas
  m_guesses<-c(alphas[1])
  for(itr in 1:iterations){
    curalphas<-E(newalphas, sources)
    newalphas<-M_basic(curalphas, sources, sink)
    m_guesses<-c(m_guesses, newalphas[1])

    if(abs(m_guesses[length(m_guesses)]-m_guesses[length(m_guesses)-1])<=10^-6) break
  }
  toret<-c(newalphas)
  return(toret)
}

"source_process_nounknown" <- function(train, envs, rarefaction_depth=1000){

  train <- as.matrix(train)

  # enforce integer data
  if(sum(as.integer(train) != as.numeric(train)) > 0){
    stop('Data must be integral. Consider using "ceiling(datatable)" or ceiling(1000*datatable) to convert floating-point data to integers.')
  }
  envs <- factor(envs)
  train.envs <- sort(unique(levels(envs)))

  # rarefy samples above maxdepth if requested
  if(!is.null(rarefaction_depth) && rarefaction_depth > 0) train <- rarefy(train, rarefaction_depth)

  # get source environment counts
  # sources is nenvs X ntaxa
  X <- t(sapply(split(data.frame(train), envs), colSums))

  rownames(X) <- c(train.envs)
  X <- t(as.matrix(X))

  return(X)
}

"read_pseudo_data"<-function(dataset){
  path_to_data<-"../data/"
  if(dataset=="DA"){
    df<-read.table(paste0(path_to_data,"DA_99_T_d10000_date_nan.txt"), fill = NA)
    return(df[complete.cases(df),])
  }else if(dataset=="DB"){
    df<-read.table(paste0(path_to_data,"DB_99_T_d10000_date_nan.txt"), fill = NA)
    return(df[complete.cases(df),])
  }else if (dataset=="F4"){
    df<-read.table(paste0(path_to_data,"F4_99_T_d10000_date_nan.txt"), fill = NA)
    return(df[complete.cases(df),])
  }else{
    df<-read.table(paste0(path_to_data,"M3_99_T_d10000_date_nan.txt"), fill = NA)
    return(df[complete.cases(df),])}
}

create_m <- function(num_sources, n, EPSILON){


  if( n == 1 ){

    index = sample(c(1:num_sources), 1)
    m_1 = runif(min = 0.6, max = 0.9, n = 1)
    resid = 1-m_1
    other_ms = resid/(num_sources-1)
    m = rep(NA, num_sources)
    m[index] = c(m_1)
    m[is.na(m)] = other_ms

  }


  if( n == 2 ){

    index = sample(c(1:num_sources), 2)
    m_1 = runif(min = 0.1, max = 0.2, n = 1)
    m_2 = runif(min = 0.4, max = 0.5, n = 1)
    resid = 1-(m_1+m_2)
    other_ms = resid/(num_sources-2)
    m = rep(NA, num_sources)
    m[index] = c(m_1, m_2)
    m[is.na(m)] = other_ms

  }


  if( n == 3 ){

    index = sample(c(1:num_sources), 3)
    m_1 = runif(min = 0.1, max = 0.5, n = 1)
    m_2 = runif(min = 0.2, max = 0.25, n = 1)
    m_3 = runif(min = 0.1, max = 0.15, n = 1)
    resid = 1-(m_1+m_2+m_3)
    other_ms = runif(min = 0.001, max = resid/(num_sources-3), n = (num_sources-3))
    m = rep(NA, num_sources)
    m[index] = c(m_1, m_2, m_3)
    m[is.na(m)] = other_ms
    m = m/sum(m)

  }
  subsum = 0
  idx = 1:length(m)

  while ((subsum+0.001) < EPSILON){
    tosub = EPSILON - subsum
    tosub = tosub / (num_sources+1)
    mask = m > tosub
    m[mask] = m[mask] - tosub
    subsum = subsum + length(m[mask]) * tosub

  }
  m = c(m,(EPSILON))

  # sum(m)
  return(m)

}


unknown_initialize <- function(sources, sink, n_sources){

  unknown_source = rep(0, length(sink))
  sum_sources = apply(sources, 2, sum)

  unknown_source = c()

  for(j in 1:length(sum_sources)){

    unknown_source[j] = max(sink[j]-sum_sources[j], 0)

  }



  return(unknown_source)

}


unknown_initialize_1 <- function(sources, sink, n_sources){

  unknown_source = rep(0, length(sink))
  sources_sum = apply(sources, 2 ,sum)


  unknown_source = c()

  for(j in 1:length(sources_sum)){

    unknown_source[j] = max(sink[j]-sources_sum[j], 0)

  }

  #Select the cor OTUs
  ind_cor = list()
  ind_known_source_abun = c()
  ind_cor_all = which(sources[1,] > 0)

  counter = matrix(0, ncol = dim(sources)[2], nrow =  dim(sources)[1])


  for(j in 1:n_sources){

    ind_cor[[j]] = which(sources[j,] > 0)

    for(k in 1:length(sources[j,])){

      if(sources[j,k] > 0){

        counter[j,k] = counter[j,k]+1
      }


    }

  }

  OTU_present_absent = apply(counter, 2, sum)
  ind_cor_all = which(OTU_present_absent >= round(n_sources*0.8))

  if(length(ind_cor_all) > 1){

    cor_abundance = round(apply(sources[,ind_cor_all], 2, median)/2) #take the min abundnace of the 'cor'
    unknown_source[ind_cor_all] = cor_abundance

  }


  #keep the sink abundance where there is no known source
  ind_no_known_source_abun = which(sources_sum == 0)

  for(j in 1:length(ind_no_known_source_abun)){

    # unknown_source[ind_no_known_source_abun[j]] = max(runif(n = 1, min = 1, max = 100), sink[ind_no_known_source_abun[j]])
    unknown_source[ind_no_known_source_abun[j]] = max((sink[ind_no_known_source_abun[j]] - rpois(n = 1, lambda = 0.5)), 0)

  }



  return(unknown_source)

}

unknown__initialize_1 <- function(sources, sink, n_sources){

  unknown_source = rep(0, length(sink))

  #zero all the OTUs with at least 1 known source
  sources_sum = apply(sources, 2 ,sum)
  ind_known_source_abun = which(sources_sum > 0)
  unknown_source[ind_known_source_abun] = 0


  #Select the cor OTUs
  ind_cor = list()
  ind_known_source_abun = c()
  ind_cor_all = which(sources[1,] > 0)

  counter = matrix(0, ncol = dim(sources)[2], nrow =  dim(sources)[1])


  for(j in 1:n_sources){

    ind_cor[[j]] = which(sources[j,] > 0)

    for(k in 1:length(sources[j,])){

      if(sources[j,k] > 0){

        counter[j,k] = counter[j,k]+1
      }


    }

  }

  OTU_present_absent = apply(counter, 2, sum)
  ind_cor_all = which(OTU_present_absent >= round(n_sources*0.8))

  if(length(ind_cor_all) > 1){

    cor_abundance = apply(sources[,ind_cor_all], 2, median) #take the median abundnace of the 'cor'
    unknown_source[ind_cor_all] = cor_abundance

  }



  #keep the sink abundance where there is no known source
  ind_no_known_source_abun = which(sources_sum == 0)

  for(j in 1:length(ind_no_known_source_abun)){

    unknown_source[ind_no_known_source_abun[j]] = max( round(sink[ind_no_known_source_abun[j]]+ rnorm(n = length(sink[ind_no_known_source_abun[j]]))), 0)

  }



  return(unknown_source)

}


FEAST <- function(source = sources_data,
                  sinks = sinks,
                  em_itr = 1000,
                  env = rownames(sources_data),
                  include_epsilon = TRUE,
                  COVERAGE,
                  unknown_initialize = 0){


  tmp = source
  test_zeros = apply(tmp, 1, sum)
  ind_to_use = as.numeric(which(test_zeros > 0))
  ind_zero = as.numeric(which(test_zeros == 0))

  source = tmp[ind_to_use,]
  sinks = sinks



  #####adding support for multiple sources#####
  totalsource<-source
  totalsource<-as.matrix(totalsource)
  sources <- split(totalsource, seq(nrow(totalsource)))
  sources<-lapply(sources, as.matrix)
  dists<-lapply(sources, function(x) x/(sum(colSums(x))))
  totaldist<-t(Reduce("cbind", dists))
  sinks<-matrix(sinks, nrow = 1, ncol = dim(totalsource)[2])

  num_sources = dim(source)[1]
  envs_simulation = c(1:(num_sources))

  source_old = source
  totalsource_old = totalsource

  source_old=lapply(source_old,t)
  source_old<- split(totalsource_old, seq(nrow(totalsource_old)))
  source_old<-lapply(source_old, as.matrix)

  #Creating the unknown source per mixing iteration
  if(include_epsilon == TRUE){

    ##Adding the initial value of the unknown source for CLS and EM
    source_2 = list()
    totalsource_2 = matrix(NA, ncol = dim(totalsource_old)[2], nrow = ( dim(totalsource_old)[1] + 1))

    for(j in 1:num_sources){

      source_2[[j]] = source_old[[j]]
      totalsource_2[j,] = totalsource_old[j,]
    }

    #create unknown for each sink i



    sinks_rarefy = rarefy(matrix(sinks, nrow = 1), maxdepth = apply(totalsource_old, 1, sum)[1]) #make

    if(unknown_initialize == 1)
      unknown_source_1 = unknown_initialize_1(sources = totalsource[c(1:num_sources),], sink = as.numeric(sinks),
                                              n_sources = num_sources)


    if(unknown_initialize == 0)
      unknown_source_1 = unknown_initialize(sources = totalsource[c(1:num_sources),], sink = as.numeric(sinks),
                                            n_sources = num_sources)

    unknown_source = unknown_source_1 + rpois(n = length(sinks), lambda = 0.5)

    unknown_source_rarefy = rarefy(matrix(unknown_source, nrow = 1), maxdepth = COVERAGE)
    source_2[[j+1]] = t(unknown_source_rarefy)
    totalsource_2[(j+1),] = t(unknown_source_rarefy)
    totalsource = totalsource_2

    source=lapply(source_2,t)
    # totalsource <- rarefy(x = totalsource, maxdepth = COVERAGE)
    source<- split(totalsource, seq(nrow(totalsource_2)))
    source<-lapply(source_2, as.matrix)

    envs_simulation <- c(1:(num_sources+1))

  }


  samps <- source
  samps<-lapply(samps, t)

  observed_samps <- samps
  observed_samps[[(num_sources + 1)]] = t(rep(0, dim(samps[[1]])[2]))


  #Calculate JSD value
  # x <- totalsource[c(1:num_sources),]
  # JSDMatrix <- jsdmatrix(x)
  # JSDMatrix <- JSDMatrix/COVERAGE
  # JS = mean(JSDMatrix[-which(JSDMatrix == 0)])
  # js_values = append(js_values, JS)
  # print(js_values)

  initalphs<-runif(num_sources+1, 0.0, 1.0)
  initalphs=initalphs/Reduce("+", initalphs)
  sink_em = as.matrix(sinks)
  pred_em<-do_EM_basic(alphas=initalphs, sources=samps, sink=sink_em, iterations=em_itr)

  tmp<-do_EM(alphas=initalphs, sources=samps, sink=sink_em, iterations=em_itr, observed=observed_samps)
  pred_emnoise = tmp$toret

  k = 1
  pred_emnoise_all = c()
  pred_em_all = c()

  for(j in 1:length(env)){

    if(j %in% ind_to_use){

      pred_emnoise_all[j] = pred_emnoise[k]
      pred_em_all[j] = pred_em[k]
      k = k+1

    }

    else{

      pred_emnoise_all[j] = 0
      pred_em_all[j] = 0
    }

  }

  pred_emnoise_all[j+1] = pred_emnoise[k]
  pred_em_all[j+1] = pred_em[k]



  names(pred_emnoise_all) = c(env,"unknown")
  names(pred_em_all) = c(env,"unknown")


  Results = list(unknown_source = unknown_source, unknown_source_rarefy = unknown_source_rarefy,
                 data_prop = data.frame(pred_emnoise_all,pred_em_all))
  return(Results)

}



#-------------
# source tracker result to plot
# The function named 'Plot_FEAST'
# You can learn more about package at:
#   https://github.com/microbiota/amplicon

#' @title FEAST Microbial traceability analysis result to plot
#' @description Input otutab, metadata or phyloseq object; ; output a table object.
#' @param result output object of the FEAST
#' @details
#' By default, input phyloseq object include metadata and otutab
#' @return  plot
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Yong-Xin Liu \email{yxliu@@genetics.ac.cn}
#' @references
#'
#' Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai.
#' NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
#' Nature Biotechnology, 2019(37), 6:676-684, DOI: \url{https://doi.org/10.1038/s41587-019-0104-4}
#'
#' @seealso beta_pcoa beta_cpcoa
#' @examples
#'
#'  Plot_FEAST(data = result)
#'
#'
#' @export
#'


# -例子
# Plot_FEAST(data = result)







# MuiPlot_FEAST(data = result)

# MuiPlot_FEAST = function(data = result){
#
#   par(mfrow=c(2,dim(result)[2]/2), mar=c(1,1,1,1))
#   # layouts = as.character(unique(metadata$SampleType))
#
#   for (i in 1:length(colnames(result))) {
#
#     labs <- paste0(row.names(result)," \n(", round(result[,i]/sum(result[,i])*100,2), "%)")
#
#     pie(result[,i],labels=labs, init.angle=90,col =  brewer.pal(nrow(result), "Reds"),
#         border="black",main =colnames(result)[i] )
#   }
# }
#



