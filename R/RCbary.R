


# result = RCbary(ps = ps ,group  = "Group",num = 3,thread = 1)
#
# RCbary = result[[1]]
# head(RCbary)
#
# path = "./Result/bNTI/"
# # dir.create(path)
# filename = paste(path,"/RCb.csv",sep = "")
# write.csv(RCbary,filename)



RCbary = function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = NULL,group  = "Group",num = 99,thread = 1){
  ps_sub <- ps
  #----------------整理map文件
  map = as.data.frame(sample_data(ps_sub))
  map$ID = row.names(map)
  sample_data(ps) = map
  #-------------------准备OTU表格
  #-----------------抽平-不设置抽平条数，默认按照最小序列数数目抽平
  set.seed(72)  # setting seed for reproducibility
  psrare = rarefy_even_depth(ps_sub )
  #检查序列数量
  sample_sums(psrare)
  # 标准化数据
  ps.norm = transform_sample_counts(psrare, function(x) x/sum(x))




  #--------------两个函数
  # 对模拟群落计算距离
  RCbray_null_func <- function(i, freq.abd.df, alpha1, alpha2, N){
    # Get simulated communities and distance
    ## initally select OTUs weighted by their frequency. The number of OTUs selected should equal the richness of the samples.
    simcom1 = data.frame(table(sample(freq.abd.df$OTU, size=alpha1, replace=F, prob=freq.abd.df$freq)), stringsAsFactors = F)
    colnames(simcom1) = c("OTU","simcom1")
    simcom1$OTU = as.character(simcom1$OTU)
    simcom1 = inner_join(simcom1, freq.abd.df, by="OTU")
    simcom2 = data.frame(table(sample(freq.abd.df$OTU, size=alpha2, replace=F, prob=freq.abd.df$freq)), stringsAsFactors = F)
    colnames(simcom2) = c("OTU","simcom2")
    simcom2$OTU = as.character(simcom2$OTU)
    simcom2 = inner_join(simcom2, freq.abd.df, by="OTU")

    ## Now recruit OTUs based on their abundance in the metacommunity
    simcom1.abd = data.frame(table(sample(simcom1$OTU, size=N-alpha1, replace=T, prob=simcom1$p)), stringsAsFactors = F)
    colnames(simcom1.abd) = c("OTU","simcom1.abd")
    simcom1.abd$OTU = as.character(simcom1.abd$OTU)
    simcom1 = full_join(simcom1, simcom1.abd, by="OTU") %>%
      mutate(simcom1.abd = ifelse(is.na(simcom1.abd), 1, simcom1.abd)) %>%
      select(OTU, simcom1.abd)

    simcom2.abd = data.frame(table(sample(simcom2$OTU, size=N-alpha2, replace=T, prob=simcom2$p)), stringsAsFactors = F)
    colnames(simcom2.abd) = c("OTU","simcom2.abd")
    simcom2.abd$OTU = as.character(simcom2.abd$OTU)
    simcom2 = full_join(simcom2, simcom2.abd, by="OTU") %>%
      mutate(simcom2.abd = ifelse(is.na(simcom2.abd), 1, simcom2.abd)) %>%
      select(OTU, simcom2.abd)


    simcom = full_join(simcom1, simcom2, by="OTU")
    simcom[is.na(simcom)] = 0
    rownames(simcom) = simcom$OTU
    simcom$OTU = NULL

    null.dist = vegdist(t(simcom), method="bray")[1]
    return(null.dist)
  }

  # 计算RCbray的主功能
  Calc_RCbray <- function(physeq, reps, nproc){
    # Get OTU table from phyloseq object
    otu.table = otu_table(physeq)

    # Get alpha diversity for each sample
    otu.PA.table = otu.table
    otu.PA.table[otu.PA.table > 0] = 1
    alpha.df = data.frame(Sample_ID = colnames(otu.PA.table), OTU.n = colSums(otu.PA.table), stringsAsFactors = F)

    # Get beta diversity matrix
    beta.table = as.matrix(vegdist(t(otu.PA.table), method="bray", diag=TRUE, upper=TRUE))

    ## Get metacommunity
    # Calculate the number of individuals in the meta community (Average read depth)
    N <- mean(apply(t(otu.table), 1, sum))

    # Calculate the average relative abundance of each taxa across communities
    p.m <- apply(t(otu.table), 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N

    # Calculate the occurrence frequency of each taxa across communities
    otu.table.bi <- 1*(t(otu.table)>0)
    freq <- apply(otu.table.bi, 2, mean)
    freq <- freq[freq != 0]

    # Combine
    freq.abd.df = data.frame(p=p, freq=freq) %>%
      tibble::rownames_to_column(var="OTU") %>%
      filter(p != 0, freq != 0) %>%
      arrange(p)

    # For each pair of samples run the RCbray analysis
    comps = combn(alpha.df$Sample_ID, m=2, simplify = F)
    RCb.df = data.frame(Site1 = character(), Site2 = character(), RCb = numeric(), stringsAsFactors = F)
    for (j in seq(1, length(comps))){
      sam = comps[[j]]
      alpha1 = alpha.df[alpha.df$Sample_ID == sam[1],]$OTU.n
      alpha2 = alpha.df[alpha.df$Sample_ID == sam[2],]$OTU.n
      # Permute "reps" many times
      rep.list = seq(1, reps)
      null.list = mclapply(rep.list, RCbray_null_func, freq.abd.df=freq.abd.df, alpha1=alpha1, alpha2=alpha2, N=N, mc.cores=nproc)

      RCb = (length(null.list[null.list > beta.table[sam[1], sam[2]]]) + (0.5*length(null.list[null.list == beta.table[sam[1], sam[2]]])))/reps
      RCb = (RCb - 0.5)*2

      RCb.df = rbind(RCb.df, data.frame(Site1=sam[1], Site2=sam[2], RCb=RCb, stringsAsFactors = F))
    }

    RCb.df
    return(RCb.df)
  }


  # 运行RCbray的计算，这个运算再5个小时左右999重复
  RCb = Calc_RCbray(psrare, num, thread)

  head(RCb)

  return(list(RCb))
}



