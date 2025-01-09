#' @title Generate Venn Diagrams and Analyze Microbial Groups in Metagenome
#' @description
#' This function creates Venn diagrams to compare microbial groups based on presence/absence thresholds,
#' and performs detailed analyses for each group. The results include Venn diagrams, bar plots, and
#' boxplots for visualizing group-specific OTU data.
#'
#' @param ps A `phyloseq` object containing microbiome data.
#' @param group A character string specifying the grouping variable in the sample metadata. Default is `"Group"`.
#' @param num An integer specifying the threshold for selecting OTUs in Venn groups. Default is `6`.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{pa}{A combined bar plot visualization of microbial compositions for each Venn segment.}
#'   \item{pb}{A combined box plot visualization for the selected microbial groups in each segment.}
#'   \item{pc}{A combined relative abundance bar plot for each Venn segment.}
#'   \item{dat.f}{A list of data frames containing bar plot data for each Venn segment.}
#'   \item{dat.f2}{A list of data frames containing statistical analysis results for each Venn segment.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Extracts OTU presence/absence data for each microbial group defined in the metadata.
#'   \item Creates a Venn diagram to identify shared and unique OTUs among groups.
#'   \item For each Venn segment, generates bar plots, performs statistical analyses, and creates boxplots to visualize group-specific differences.
#'   \item Outputs the processed data and visualizations for further analysis.
#' }
#'
#' Statistical tests are performed to identify significant differences in microbiome abundances between groups for each Venn segment. The results are returned as data frames and visualized using boxplots.
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' result = VenSuper.metm(ps = ps.16s,group =  "Group",num = 6)
#'  p7_1 <- result[[1]]
#'  p7_1+ scale_fill_manual(values = colset1)+ scale_color_manual(values = colset1,guide = F)
#'  p7_2 <- result[[3]]
#'  p7_2
#'  p8 <- result[[2]]


#清空内存
# rm(list=ls())
# ps = readRDS("./ps_liu.rds")
# result =VenSeper(ps,num = 6,path = "./phyloseq_3-4_ven")
# # 提取韦恩图中全部部分的otu极其丰度做门类柱状图
# result[[1]]
# #每个部分序列的数量占比，并作差异
# result[[2]]
# # 每部分的otu门类冲积图
# result[[3]]
# #
# otu = NULL
# tax = NULL
# map = NULL
# tree = NULL
# ps = ps
# ps
# num = 6
# group = "Group"
# path = "./"



VenSuper.metm =  function(ps = NULL,group  = "Group",num = 6){


  # library (VennDiagram)
  rep = num
  ps
  mapping = as.data.frame(sample_data(ps))
  ps1 = ps
  ps1

  ps1_rela  = transform_sample_counts(ps1, function(x) x / sum(x) );ps1_rela
  print("1")


  aa = vegan_otu(ps1)
  otu_table = as.data.frame(t(aa))
  count = aa
  countA = count
  dim(count)
  sub_design <- as.data.frame(sample_data(ps1))

  sub_design $Group
  levels(sub_design $Group)

  pick_val_num <- num*2/3
  count[count > 0] <- 1

  count2 = as.data.frame(count)


  iris.split <- split(count2,as.factor(sub_design$Group))
  iris.apply <- lapply(iris.split,function(x)colSums(x[]))
  iris.combine <- do.call(rbind,iris.apply)
  ven2 = t(iris.combine)

  ven2[ven2 < pick_val_num]  = 0
  ven2[ven2 >=pick_val_num]  = 1
  ven2 = as.data.frame(ven2)


  ven3 = as.list(ven2)
  ven2 = as.data.frame(ven2)
  head(ven2)


  ven3 = as.list(ven2)
  for (i in 1:ncol(ven2)){

    ven3[[i]] <-  row.names(ven2[ven2[i] == 1,])

  }
  ven_pick = VennDiagram::get.venn.partitions(ven3)
  ven_pick = as.matrix(ven_pick)
  ven_pick = as.data.frame(ven_pick)
  # colnames(ven_pick)
  # row.names(ven_pick)


  sub_design <- as.data.frame(sample_data(ps1_rela ))


  mapping = as.data.frame(sample_data(ps))
  length(ven_pick$..set..)
  plot1 = list()
  plot2 = list()
  plot3 = list()
  dat.f = list()
  dat.f2 = list()
  #i=1
  for (i in 1:length(ven_pick$..set..) ) {

    # 查看i是韦恩图哪一个部分的otu
    ven_pick$..set..[[i]]


    abs = length(unique(mapping$Group))
    abs = 2+abs
    aab = ven_pick[[i,abs]]


    if (length(aab) == 0) {
      i = i +1
    }




    # #提取otu子集
    otu = as.data.frame(t(vegan_otu(ps1_rela)))
    dim(otu)
    otu$ID = row.names(otu)
    otu1<- dplyr::filter(otu, ID %in% aab)
    head(otu1)
    row.names(otu1) = otu1$ID

    otu1$ID = NULL
    subtab = as.matrix(otu1)
    ps_sub = ps1_rela
    print("2")
    ps_sub <- phyloseq::phyloseq(otu_table(subtab, taxa_are_rows=TRUE),
                                 phyloseq::tax_table(ps_sub),
                                 sample_data(ps_sub)
    )
    ps_sub



    #--------------过滤掉ven中没有otu的分组
    num = ven_pick[i,1:length(unique(mapping$Group))]
    num = as.data.frame(num )
    colnames(num)[num[1,] == FALSE]
    # ps_sub  <- subset_samples( ps_sub,!Group %in% colnames(num)[num[1,] == FALSE]);ps_sub
    map = as.data.frame(sample_data(ps_sub))

    map$ID = row.names(map)
    head(map)
    maps<- dplyr::filter(as.tibble(map),!Group %in% colnames(num)[num[1,] == FALSE]) %>% as.data.frame()
    row.names(maps) = maps$ID
    ps_sub = ps_sub
    sample_data( ps_sub ) = maps


    # 这部分进行堆叠柱状图的出图，但是不能再做标准化了

    result = barMainplot.micro(ps = ps_sub,j = "Phylum",axis_ord = NULL,label = FALSE ,sd = FALSE,
                               Top = 10,tran = FALSE)
    print("3")
    #提取脱图片
    p = result[[1]]
    # 提取作图数据，也就是这部分ven包含otu的数量
    result[[2]]
    #提取冲击图
    p3 = result[[3]]

    # filename = paste(path,"/TAX_ven_pick_",ven_pick$..set..[[i]],".csv",sep = "")
    # write.csv(result[[2]],filename,quote = FALSE)
    dat.f[[i]] = result[[2]]
    names(dat.f)[i] = paste("TAX_ven_pick_",ven_pick$..set..[[i]],sep = "")

    #  对ven每个部分共有或者特有的部分做差异分析，这里使用我开发了
    vencount = as.data.frame(sample_sums(ps_sub ))
    end = merge(vencount ,sub_design ,by = "row.names",all = FALSE)

    # 准备对这部分序列进行差异分析和出图
    data_wt = data.frame(ID= end$Row.names,group = end$Group,count = end$`sample_sums(ps_sub)`)
    head(data_wt)
    colnames(data_wt)[3] =ven_pick$..set..[[i]]
    # library(EasyAovWlxPlot)
    if (length(unique(data_wt$group)) !=1) {
      result = KwWlx2(data = data_wt, i= 3)
      PlotresultBox = aovMuiBoxP2(data = data_wt, i= 3,sig_show ="abc",result = result[[1]])
      # 提取检验结果
      PlotresultBox[[2]]
      #提取图片
      p2 = PlotresultBox[[1]]
      p2

    }

    if (length(unique(data_wt$group)) ==1) {
      Mytheme <- theme_bw()+

        # scale_fill_manual(values = mi, guide = guide_legend(title = NULL))+
        theme(

          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),

          plot.title = element_text(vjust = -8.5,hjust = 0.1),
          axis.title.y =element_text(size = 20,face = "bold",colour = "black"),
          axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
          axis.text = element_text(size = 20,face = "bold"),
          axis.text.x = element_text(colour = "black",size = 14),
          axis.text.y = element_text(colour = "black",size = 14),
          legend.text = element_text(size = 15,face = "bold"),
          legend.position = "none"#是否删除图例

        )
      aa  =colnames(data_wt)[3]
      colnames(data_wt)[3] = "XX"

      p2 = ggplot(data_wt) + geom_boxplot(aes(x = group, y =XX,color = group))  + Mytheme  +
        labs(x="",
             y=aa)

    }



    # filename = paste(path,"/SeqStat_ven_pick_",ven_pick$..set..[[i]],".csv",sep = "")
    # write.csv(PlotresultBox[[2]],filename,quote = FALSE)
    dat.f2[[i]] = PlotresultBox[[2]]
    names(dat.f2)[i]  = paste("SeqStat_ven_pick_",ven_pick$..set..[[i]],sep = "")

    plot1[[i]] =p
    plot2[[i]] =p2
    plot3[[i]] =p3


  }


  pa  = ggpubr::ggarrange(plotlist = plot1, common.legend = TRUE, legend="right")
  pa


  pb  = ggpubr::ggarrange(plotlist = plot2, common.legend = TRUE, legend="right")
  pb


  pc  = ggpubr::ggarrange(plotlist = plot3, common.legend = TRUE, legend="right")
  pc


  return(list(pa,pb,pc, dat.f,dat.f2))

}


