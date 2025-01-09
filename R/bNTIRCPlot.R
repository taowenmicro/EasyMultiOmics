
#' @title Integrating bNTI and RCbray for Ecological Processes Visualization
#'
#' @description
#' The `bNTIRCPlot` function combines β-nearest taxon index (bNTI) and Raup-Crick dissimilarity with Bray-Curtis distance (RCbray) to visualize the relative contributions of ecological processes shaping microbial communities.
#'
#' @param otu A data frame or matrix containing OTU (Operational Taxonomic Unit) abundance data. Rows represent taxa, and columns represent samples.
#' @param tax A data frame containing taxonomic annotations for each OTU.
#' @param map A data frame containing sample metadata, including group information.
#' @param tree A phylogenetic tree object.
#' @param ps A `phyloseq` object containing OTU, taxonomic, and sample data. If provided, this supersedes `otu`, `tax`, `map`, and `tree`.
#' @param RCb A data frame containing RCbray values for each pair of samples.
#' @param bNTI A data frame containing bNTI values for each pair of samples.
#' @param group A character string specifying the grouping variable in the metadata. Default is `"Group"`.
#'
#' @return A list containing:
#' \itemize{
#'   \item  A boxplot of bNTI values within each group.
#'   \item  A bar plot of the relative proportions of ecological processes for each group.
#'   \item `A combined plot of the bNTI and ecological process bar plot.
#'   \item  A merged data frame containing bNTI, RCbray, and ecological process information.
#'   \item  A summary table of the proportions of ecological processes for each group.
#' }
#'
#' @details
#' This function integrates bNTI and RCbray metrics to classify ecological processes into five categories:
#' \itemize{
#'   \item `Drift`
#'   \item `Dispersal Limited`
#'   \item `Homogenizing Dispersal`
#'   \item `Variable Selection`
#'   \item `Homogeneous Selection`
#' }
#' It generates visualizations to depict the contribution of these processes across microbial communities.
#'
#' @examples
#' \dontrun{
#' result = bNTICul(ps = psphy,group  = "Group",num = 10,thread = 1)
#' bNTI = result[[1]]
#' result = RCbary(ps = psphy ,group  = "Group",num = 10,thread = 1)
#' RCbary = result[[1]]
#' result = bNTIRCPlot(ps = psphy ,RCb  =Rcb, bNTI = bNTI,group  = "Group")
#' }
#'
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @export
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
