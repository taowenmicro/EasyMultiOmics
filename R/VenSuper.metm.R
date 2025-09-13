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
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njau.edu.cn}
#' @examples
#' result = VenSuper.metm(ps = ps.16s,group =  "Group",num = 6)
#'  p7_1 <- result[[1]]
#'  p7_1+ scale_fill_manual(values = colset1)+ scale_color_manual(values = colset1,guide = F)
#'  p7_2 <- result[[3]]
#'  p7_2
#'  p8 <- result[[2]]
#' @export


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


VenSuper.metm <- function(ps, group = "Group", num = 6) {

  # Step1: 基础数据准备
  ps_rela <- phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
  otu_tab <- as.data.frame(t(ggClusterNet::vegan_otu(ps_rela)))
  mapping <- as.data.frame(phyloseq::sample_data(ps_rela))

  # Step2: OTU presence/absence
  count <- ggClusterNet::vegan_otu(ps)
  count[count > 0] <- 1
  count_df <- as.data.frame(count)

  # 按组汇总OTU
  iris_split <- split(count_df, as.factor(mapping[[group]]))
  iris_apply <- lapply(iris_split, function(x) colSums(x))
  iris_combine <- do.call(rbind, iris_apply)
  ven2 <- t(iris_combine)

  pick_val_num <- num * 2 / 3
  ven2[ven2 < pick_val_num] <- 0
  ven2[ven2 >= pick_val_num] <- 1
  ven2 <- as.data.frame(ven2)

  # Step3: 生成Venn partitions
  ven_list <- lapply(1:ncol(ven2), function(i) rownames(ven2[ven2[, i] == 1, ]))
  names(ven_list) <- colnames(ven2)
  ven_pick <- VennDiagram::get.venn.partitions(ven_list)
  ven_pick <- as.data.frame(as.matrix(ven_pick))

  # Step4: 初始化输出列表
  plot_bar_list <- list()
  plot_box_list <- list()
  plot_stack_list <- list()
  dat_f_list <- list()
  dat_f2_list <- list()

  # Step5: 循环处理每个 Venn 分区
  # i = 1
  for(i in seq_len(nrow(ven_pick))) {

    # 当前Venn区 OTU
    aab <- unlist(ven_pick[i, (2 + length(unique(mapping[[group]]))):ncol(ven_pick)])
    aab <- aab[aab != ""]
    if(length(aab) == 0) next

    # 构建子集phyloseq对象
    otu_sub <- otu_tab[rownames(otu_tab) %in% aab, , drop = FALSE]
    ps_sub <- phyloseq::phyloseq(
      phyloseq::otu_table(as.matrix(otu_sub), taxa_are_rows = TRUE),
      phyloseq::tax_table(ps_rela),
      phyloseq::sample_data(ps_rela)
    )

    # 去掉没有OTU的分组
    sub_groups <- colnames(ven2)[which(ven2[i, ] != 0)]
    sample_map <- as.data.frame(phyloseq::sample_data(ps_sub))
    sample_map <- sample_map[sample_map[[group]] %in% sub_groups, , drop = FALSE]
    rownames(sample_map) <- rownames(sample_map)
    phyloseq::sample_data(ps_sub) <- sample_map

    # Step6: 堆叠柱状图
    bar_res <- barMainplot.metm(ps = ps_sub, j = "Phylum", Top = 10, tran = FALSE, sd = FALSE, label = FALSE)
    plot_bar_list[[i]] <- bar_res[[1]]
    plot_stack_list[[i]] <- bar_res[[3]]
    dat_f_list[[i]] <- bar_res[[2]]
    names(dat_f_list)[i] <- paste0("TAX_ven_pick_", ven_pick$..set..[i])

    # Step7: 差异分析与boxplot
    vencount <- as.data.frame(sample_sums(ps_sub))
    colnames(vencount) <- "count"
    # vencount$group <- as.character(mapping[rownames(vencount), group])
    idx <- match(rownames(vencount), rownames(mapping))
    group_labels <- as.character(mapping[[group]][idx])
    vencount$group <- group_labels
    if(length(unique(vencount$group)) > 1) {
      kw_res <- KwWlx2(data = data.frame(ID = rownames(vencount), group = vencount$group, count = vencount$count), i = 3)
      box_res <- aovMuiBoxP2(data = data.frame(ID = rownames(vencount), group = vencount$group, count = vencount$count),
                             i = 3, sig_show = "abc", result = kw_res[[1]])
      plot_box_list[[i]] <- box_res[[1]]
      dat_f2_list[[i]] <- box_res[[2]]
    } else {
      # 单组情况直接绘制boxplot
      df_single <- data.frame(group = vencount$group, count = vencount$count)
      plot_box_list[[i]] <- ggplot2::ggplot(df_single, ggplot2::aes(x = group, y = count, color = group)) +
        ggplot2::geom_boxplot() +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "", y = ven_pick$..set..[i])
      dat_f2_list[[i]] <- data.frame()
    }
    names(dat_f2_list)[i] <- paste0("SeqStat_ven_pick_", ven_pick$..set..[i])
  }

  # Step8: 组合输出
  pa <- ggpubr::ggarrange(plotlist = plot_bar_list, common.legend = TRUE, legend = "right")
  pb <- ggpubr::ggarrange(plotlist = plot_box_list, common.legend = TRUE, legend = "right")
  pc <- ggpubr::ggarrange(plotlist = plot_stack_list, common.legend = TRUE, legend = "right")

  return(list(pa, pb, pc, dat_f_list, dat_f2_list))
}

