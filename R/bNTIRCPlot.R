

#---
#
# # bNTI data
# # Import bNTI data
#
# library(EasyMicrobiome)
# library(dplyr)
# library(ggplot2)
# bNTI = read.csv("../pipeline//bNTI.csv",row.names = 1)
# head(bNTI)
#
# ps = readRDS("../ori_data/ps_liu.rds")
#
# # RCbray 数据读入，修改列名
# RCb = read.csv("../pipeline/RCb.csv",row.names = 1) %>%
#   mutate(Sample_1 = Site2, Sample_2 = Site1)
# head(RCb)
# result = bNTIRCPlot(ps = ps ,RCb  = RCb,bNTI = bNTI,group  = "Group")
#
# #--bNTI出图片
# result[[1]]
#
# #RCbary可视化
# result[[2]]
# #组合图片BNTI，RCbray
# result[[3]]
#
# plotdata = result[[4]]
# head(plotdata)


bNTIRCPlot = function(otu = NULL,tax = NULL,
                      map = NULL,tree = NULL ,
                      ps = NULL,
                      RCb  = RCb,bNTI = bNTI,group  = "Group"){

  ps = inputMicro(otu,tax,map,tree,ps,group  = group)
  ps

  psrare <- ps
  map = as.data.frame(sample_data(psrare))
  map$ID = row.names(map)
  sample_data(psrare) = map

  # Get habitat metadata and add it to the βNTI then merge with the RCbray dataset
  eco.meta1 = data.frame(sample_data(psrare)) %>%
    select(ID, Group) %>%
    dplyr::rename(Sample_1 = ID, Group_1 = Group)

  eco.meta2=data.frame(sample_data(psrare)) %>%
    select(ID, Group) %>%
    dplyr::rename(Sample_2 = ID, Group_2 = Group)

  # bNTI 匹配第一列和第二列的分组信息
    bNTI.df = inner_join(bNTI, eco.meta1) %>%
    inner_join(eco.meta2)

  # 合并两个数据
  turnover.df = inner_join(bNTI.df, RCb)
  head(turnover.df)
  dim(turnover.df)


  #--------------合并文件保存
  # write.csv(turnover.df,"./Result/bNTI//bNTI_RCbray.csv")



  #-----按照分组统计作图

  #------------bNIT作图
  dim(bNTI.df)
  within.bNTI.df = bNTI.df %>%
    filter(Group_1 == Group_2) %>%
    mutate(Group = Group_1)

  head(within.bNTI.df )

  # map$Group
  # ecosystem.conv = data.frame(Group = c("KO", "OE", "WT"), Group2 = c("Cropland", "Old-field", "Forest"))
  # within.bNTI.df = left_join(within.bNTI.df, ecosystem.conv)
  # within.bNTI.df$Group2 = factor(within.bNTI.df$Group2, levels = c("Cropland", "Old-field", "Forest"))
  # within.bNTI.df$Group = factor(within.bNTI.df$Group, levels=c("KO", "OE", "WT"))
  # eco.bNTI.plot = ggplot(within.bNTI.df, aes(x=Group, y=bNTI,fill = "Group")) +
  #   geom_boxplot(outlier.shape=1) +
  #   geom_hline(yintercept = 2, linetype=2, size=0.5) +
  #   geom_hline(yintercept = -2, linetype=2, size=0.5) +
  #   labs(x="", y="bNTI") +
  #   theme_bw() +
  #   theme(legend.position = "none",
  #         axis.text = element_text(size=12),
  #         axis.text.x = element_text(angle=45, hjust=1),
  #         axis.title = element_text(size=14))
  # eco.bNTI.plot
  eco.bNTI.plot <- ggplot(within.bNTI.df, aes(x=Group, y=bNTI)) +
    geom_jitter(alpha = 0.1,color ="#984EA3") +
    geom_boxplot(outlier.shape=1,outlier.alpha = 0,fill = "#984EA3") +

    geom_hline(yintercept = 2, linetype=2, size=0.5) +
    geom_hline(yintercept = -2, linetype=2, size=0.5) +
    labs(x="", y="bNTI") +
    theme_classic() +
    theme(legend.position = "none",
          axis.text = element_text(size=12),
          axis.text.x = element_text(angle=45, hjust=1),
          axis.title = element_text(size=14))



  # 现在按照RCbray进行分开标记系统发育过程
  eco.turnover.df = turnover.df %>%
    filter(Group_1 == Group_2) %>%
    mutate(Group = Group_1)

  # eco.turnover.df = left_join(eco.turnover.df, ecosystem.conv)
  # eco.turnover.df$Group2 = factor(eco.turnover.df$Group2, levels = c("Cropland", "Old-field", "Forest"))
  # eco.turnover.df$Group = factor(eco.turnover.df$Group, levels=c("KO", "OE", "WT"))


  head(eco.turnover.df )


  ## Calculate the relative influence of each process
  eco.turnover.df = eco.turnover.df %>%
    mutate(process = ifelse(abs(bNTI) < 2,
                            ifelse(abs(RCb) < 0.95, "Drift",
                                   ifelse(RCb >= 0.95, "Dispersal Limited",
                                          ifelse(RCb <= -0.95, "Homogenizing Dispersal", "ERROR"))),
                            ifelse(bNTI >= 2, "Variable Selection",
                                   ifelse(bNTI <= -2, "Homogeneous Selection", "ERROR"))))


  eco.turnover.df$process = factor(eco.turnover.df$process, levels = c("Drift",
                                                                       "Dispersal Limited", "Homogenizing Dispersal",
                                                                       "Variable Selection", "Homogeneous Selection"))

  head(eco.turnover.df)


  #------计算每个组的系统发育过程中五个部分分别占有的比例
  pre = eco.turnover.df %>%
    dplyr::group_by(Group, process) %>%
    dplyr::summarize(n_sites = n(),
                     perc=(n()/45)*100) %>%
    as.data.frame
  # head(numeco  )
   numeco <- pre %>%  dplyr::group_by(Group) %>%
     dplyr::summarise(num = sum(n_sites))
   alleco <- pre %>% dplyr::left_join(numeco,by = "Group")
   alleco$perc =  alleco$n_sites/ alleco$num * 100
   sum.eco.turnover.df = alleco
  eco.turnover.plot = ggplot(sum.eco.turnover.df, aes(x=Group, y=perc, fill=process)) +
    geom_bar(stat="identity", color="black") +
    # scale_fill_manual(values = c("white", "grey75", "grey50", "black")) +
    labs(x="", y="Percent of pairs (%)", fill="Process") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size=12),
          axis.text.x = element_text(angle=45, hjust=1),
          axis.title = element_text(size=14),
          legend.key.size = unit(10, "mm"),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14))
  eco.turnover.plot


  # Merge the plots
  eco.plot = cowplot::plot_grid(eco.bNTI.plot, eco.turnover.plot,
                                rel_widths=c(0.6, 1), labels=c("A", "B"))
  eco.plot


  return(list( eco.bNTI.plot, eco.turnover.plot,eco.plot,turnover.df,sum.eco.turnover.df))
}
