#' @title Generate Circular Plot for Microbial Composition
#'
#' @description
#' This function creates a circular plot to visualize microbial community composition based on OTU abundance.
#' It summarizes the microbial taxa at the specified rank, calculates the mean abundance within groups, and displays the data in a circular chord diagram.
#'
#' @param ps A `phyloseq` object containing microbiome data.
#' @param Top An integer specifying the number of top taxa to display. Taxa not in the top are grouped as `"others"`. Default is `10`.
#' @param rank A numeric or character value specifying the taxonomic rank for aggregation (e.g., `"Phylum"`, `"Genus"`, etc.). Default is `7`.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{mer_otu_mean}{A matrix summarizing the mean abundance of OTUs across groups.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Transforms the data into relative abundances using `phyloseq::transform_sample_counts`.
#'   \item Aggregates microbial taxa at the specified rank using `ggClusterNet::tax_glom_wt`.
#'   \item Calculates the mean abundance of OTUs for each group and selects the top `Top` taxa.
#'   \item Groups less abundant taxa into a category labeled `"others"`.
#'   \item Visualizes the relationships between microbial taxa and groups using a circular chord diagram with the `circlize` package.
#' }
#'
#' @examples
#' \dontrun{
#' res = cir_plot.micro(ps  = ps.16s,Top = 12,rank = 6)
#'
#' }
#'
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#'
#' @export


cir_plot.micro = function(ps  =ps,
                    Top = 10,
                    rank = 7
){
  ps_rela = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) );ps_rela

  ps_P <- ps_rela %>%
    ggClusterNet::tax_glom_wt( rank = rank)
  ps_P

  otu_P = as.data.frame((ggClusterNet::vegan_otu(ps_P)))
  head(otu_P)
  tax_P = as.data.frame(ggClusterNet::vegan_tax(ps_P))

  sub_design <- as.data.frame(phyloseq::sample_data(ps_P))
  count2 =   otu_P
  #数据分组
  iris.split <- split(count2,as.factor(sub_design$Group))
  #数据分组计算平均值
  iris.apply <- lapply(iris.split,function(x)colSums(x[]))
  # 组合结果
  iris.combine <- do.call(rbind,iris.apply)
  ven2 = t(iris.combine)
  # head(ven2)

  if (is.numeric(rank)) {
    lev = phyloseq::rank_names(ps)[rank]
  } else{
    lev = rank
  }


  Taxonomies <- ps %>%
    ggClusterNet::tax_glom_wt(rank = rank) %>%
    phyloseq::transform_sample_counts(function(x) {x/sum(x)} )%>%
    phyloseq::psmelt() %>%
    #dplyr::filter(Abundance > 0.05) %>%
    dplyr::arrange( !!sym(lev))
  iris_groups<- dplyr::group_by(Taxonomies, !!sym(lev))
  ps0_sum <- dplyr::summarise(iris_groups, mean(Abundance), sd(Abundance))
  ps0_sum[is.na(ps0_sum)] <- 0
  head(ps0_sum)
  colnames(ps0_sum) = c("ID","mean","sd")


  ps0_sum <- dplyr::arrange(ps0_sum,desc(mean))
  ps0_sum$mean <- ps0_sum$mean *100
  ps0_sum <- as.data.frame(ps0_sum)
  head(ps0_sum)
  top_P = ps0_sum$ID[1:Top];top_P

  ### 开始进一步合并过滤
  otu_P = as.data.frame(t(otu_P))
  otu_tax = merge(ven2,tax_P,by = "row.names",all = F)
  dim(otu_tax)
  otu_tax[,lev] = as.character(otu_tax[,lev])
  otu_tax[,lev][is.na(otu_tax[,lev])] = "others"

  i = 1
  for (i in 1:nrow(otu_tax)) {
    if(otu_tax[,lev] [i] %in% top_P){otu_tax[,lev] [i] = otu_tax[,lev] [i]}

    else if(!otu_tax[,lev] [i] %in% top_P){otu_tax[,lev] [i] = "others"}

  }

  otu_tax[,lev] = as.factor(otu_tax[,lev])
  head(otu_tax)

  otu_mean = otu_tax[as.character(unique(sub_design$Group))]
  head(otu_mean)
  row.names(otu_mean) = row.names(otu_tax)
  iris.split <- split(otu_mean,as.factor(otu_tax[,lev]))
  #数据分组计算平均值
  iris.apply <- lapply(iris.split,function(x)colSums(x[]))
  # 组合结果
  iris.combine <- do.call(rbind,iris.apply)
  mer_otu_mean = t(iris.combine)

  head(mer_otu_mean )


  # mer_otu_mean = t(mer_otu_mean)


  # library(statnet)
  # library(circlize)
  # library(RColorBrewer)#调色板调用包
  mi_sam = RColorBrewer::brewer.pal(9,"Set1")
  mi_tax = colorRampPalette(RColorBrewer::brewer.pal(9,"Set3"))(length(row.names(mer_otu_mean)))
  # library("scales")
  # show_col(mi_sam )
  # show_col(mi_tax)
  # circlize::CELL_META

  grid.col = NULL
  #这里设置样品颜色
  grid.col[as.character(unique(sub_design$Group))] = mi_sam
  #设置群落中物种水平颜色
  # grid.col[colnames(mer_otu_mean)] = mi_tax
  grid.col[row.names(mer_otu_mean)] = mi_tax

  #
  # FileName2 <- paste(path,"/",lev,"_cricle",".pdf", sep = "")
  # pdf(FileName2, width = 12, height = 8)

  #gap.degree修改间隔，不同小块之间的间隔
  circlize::circos.par(gap.degree = c(rep(2, nrow(mer_otu_mean)-1), 10, rep(2, ncol(mer_otu_mean)-1), 10),
             start.degree = 180)
   circlize::chordDiagram(mer_otu_mean,
               directional = F,
               diffHeight = 0.06,
               grid.col = grid.col,
               reduce = 0,
               transparency = 0.5,
               annotationTrack =c("grid", "axis"),
               preAllocateTracks = 2
  )

  circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
  circlize::circos.text(circlize::CELL_META$xcenter, circlize::CELL_META$ylim[1], circlize::CELL_META$sector.index,
  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))}, bg.border = NA)# here set bg.border to NA is important


 circlize::circos.clear()
  #dev.off()

#   grid.col = NULL
#   #这里设置样品颜色
#   grid.col[as.character(unique(sub_design$Group))] = mi_sam
#   #设置群落中物种水平颜色
#   # grid.col[colnames(mer_otu_mean)] = mi_tax
#   grid.col[row.names(mer_otu_mean)] = mi_tax
# #
# #     FileName2 <- paste(path,"/",lev,"_cricle",".png", sep = "")
# #
# #
# #   png(FileName2, res=150, width = 1000, height = 1008)
#
#   #gap.degree修改间隔，不同小块之间的间隔
#   circlize::circos.par(gap.degree = c(rep(2, nrow(mer_otu_mean)-1), 10, rep(2, ncol(mer_otu_mean)-1), 10),
#              start.degree = 180)
#   circlize::chordDiagram(mer_otu_mean,directional = F,
#                reduce = 0,
#                diffHeight = 0.06,grid.col = grid.col, transparency = 0.5, annotationTrack =c("grid", "axis"),
#                preAllocateTracks = 2
#   )
#
#   circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
#     circlize::circos.text(circlize::CELL_META$xcenter, circlize::CELL_META$ylim[1], circlize::CELL_META$sector.index,
#                 facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))}, bg.border = NA) # here set bg.border to NA is important
#   circlize::circos.clear()
#   dev.off()
  return(list(mer_otu_mean ))
}
