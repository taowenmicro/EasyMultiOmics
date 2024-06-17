# library(EasyMicrobiome)
# ps = readRDS("./ori_data/ps_liu.rds")
# ps
#
# result = neutralModel(ps = ps,group  = "Group")
# #--合并图表
# p =  result[[1]]
#
# #--分开的图表--存储在一个list中
# plist = result[[2]]
# #-提取单个的图表
# plist[[1]]

neutralModel = function(otu = NULL,
                        tax = NULL,
                        map = NULL,
                        tree = NULL,
                        ps = NULL,
                        group  = "Group",
                        ncol = 3,
                        nrow  = 1
                        
                        ){

  # 抽平，默认使用最小序列抽平
  ps = inputMicro(otu,tax,map,tree,ps,group  = group)
  ps
  set.seed(72)  #设置随机种子，保证结果可重复
  psrare = rarefy_even_depth(ps)

  # 标准化
  ps.norm = transform_sample_counts(psrare, function(x) x/sum(x))



  #------------------------------------------开始计算中性模型----------------------------------------------------------
  map = as.data.frame(sample_data(psrare))
  aa = levels(map$Group)
  aa
  map$ID = row.names(map)

  plots = list()
  dat1 = list()
  dat2 = list()
  i =1
  for (i in 1:length(aa)) {


    maps<- dplyr::filter(as.tibble(map),Group %in%aa[i])
    maps = as.data.frame(maps)
    row.names(maps) = maps$ID
    ps_sub = psrare
    sample_data( ps_sub ) =maps ;ps_sub

    # 提取OTU表格
    OTU.table = t(otu_table(ps_sub))
    head(OTU.table )
    # 将整个群落看做一个整体，计算每个样本的序列数，并求取均值Calculate the number of individuals in the meta community (Average read depth)
    N <- mean(apply(OTU.table, 1, sum))

    #计算每个OTU的的平均序列数 Calculate the average relative abundance of each taxa across communities
    p.m <- apply(OTU.table, 2, mean)
    #去除OTU序列数为0的OTU
    p.m <- p.m[p.m != 0]
    p <- p.m/N
    p.df = data.frame(p) %>%
      rownames_to_column(var="OTU")

    # Calculate the occurrence frequency of each taxa
    OTU.table.bi <- 1*(OTU.table>0)
    freq.table <- apply(OTU.table.bi, 2, mean)
    freq.table <- freq.table[freq.table != 0]
    freq.df = data.frame(OTU=names(freq.table), freq=freq.table)

    #Combine
    C <- inner_join(p.df,freq.df, by="OTU") %>%
      arrange(p)
    # Remove rows with any zero (absent in either source pool or local communities). You already did this, but just to make sure we will do it again.
    C.no0 <- C %>%
      filter(freq != 0, p != 0)

    #Calculate the limit of detection
    d <- 1/N

    ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
    p.list <- C.no0$p
    freq.list <- C.no0$freq
    m.fit <- nlsLM(freq.list ~ pbeta(d, N*m*p.list, N*m*(1-p.list), lower.tail=FALSE), start=list(m=0.1))
    m.ci <- confint(m.fit, 'm', level=0.95)
    m.sum <- summary(m.fit)
    m.coef = coef(m.fit)

    freq.pred <- pbeta(d, N*coef(m.fit)*p.list, N*coef(m.fit)*(1-p.list), lower.tail=FALSE)
    Rsqr <- 1 - (sum((freq.list - freq.pred)^2))/(sum((freq.list - mean(freq.list))^2))

    # Get table of model fit stats
    fitstats <- data.frame(m=m.coef, m.low.ci=m.ci[1], m.up.ci=m.ci[2],
                           Rsqr=Rsqr, p.value=m.sum$parameters[4], N=N,
                           Samples=nrow(OTU.table), Richness=length(p.list),
                           Detect=d)

    # Get confidence interval for predictions
    freq.pred.ci <- binconf(freq.pred*nrow(OTU.table), nrow(OTU.table), alpha=0.05, method="wilson", return.df=TRUE)

    # Get table of predictions
    pred.df <- data.frame(metacomm_RA=p.list, frequency=freq.pred,
                          frequency_lowerCI=freq.pred.ci[,2],
                          frequency_upperCI=freq.pred.ci[,3]) %>%
      unique()

    # Get table of observed occupancy and abundance
    obs.df = C.no0 %>%
      dplyr::rename(metacomm_RA = p, frequency=freq)

    head(obs.df)



    p = ggplot(data=obs.df) +
      geom_point(data=obs.df, aes(x=log10(metacomm_RA), y=frequency),
                 alpha=.3, size=2, color="#8DD3C7") +
      geom_line(data=pred.df, aes(x=log10(metacomm_RA), y=frequency), color="#FFFFB3") +
      geom_line(data=pred.df, aes(x=log10(metacomm_RA), y=frequency_lowerCI), linetype=2, color="#FFFFB3") +
      geom_line(data=pred.df, aes(x=log10(metacomm_RA), y=frequency_upperCI), linetype=2, color="#FFFFB3") +
      # geom_text(data=fitstats, aes(label = paste("R^2 == ", round(Rsqr, 3))),
      #           x=1, y=0.75, size=4, parse=TRUE) +
      # geom_text(data=fitstats, aes(label = paste("italic(m) ==", round(m, 3))),
      #           x=-1, y=0.85, size=4, parse=TRUE) +
      labs(x="Log10 abundance in\nmetacommunity", y="Frequency detected",title = paste(aa[i],paste("R^2 == ", round(fitstats$Rsqr, 3)),paste("italic(m) ==", round(fitstats$m, 3)))) +
      theme_bw() +
      theme(axis.line = element_line(color="black"),
            legend.position = "none",
            axis.title = element_text(size=14),
            axis.text = element_text(size=12))

    p



    plots[[aa[i]]] = p
    dat1[[aa[i]]] = obs.df
    dat2[[aa[i]]] = pred.df
  }


  # plots$ABCD
  # library(ggpubr)
  # nrow=2,,ncol=4
  p  = ggpubr::ggarrange(plotlist = plots,common.legend = TRUE, legend="right",ncol = ncol,nrow = nrow)
  p

  return(list(p,plots,dat1,dat2))

}
