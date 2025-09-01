#' @title Correlation Analysis and Visualization of Environmental Variables and OTUs
#'
#' @description
#' The `cor_env_ggcorplot_top` function calculates the correlation between environmental variables and top OTUs (or taxa),
#' and visualizes the results as a heatmap and scatterplot. It also supports hierarchical clustering for both rows and columns.
#'
#' @param ps A `phyloseq` object containing OTU/taxa abundance data and sample metadata.
#' @param jj A taxonomic rank (e.g., "Genus", "Family") for aggregating OTUs. If `NULL`, no aggregation is performed. Default is `NULL`.
#' @param tran A logical value indicating whether to transform OTU counts into relative abundances. Default is `TRUE`.
#' @param Top An integer specifying the number of top OTUs to retain based on total abundance. Default is `10`.
#' @param env1 A data frame of environmental variables, where rows are samples and columns are variables.
#' @param label A logical value indicating whether to label the rows and columns in the plot. Default is `TRUE`.
#' @param col_cluster A logical value indicating whether to perform hierarchical clustering on the columns. Default is `TRUE`.
#' @param row_cluster A logical value indicating whether to perform hierarchical clustering on the rows. Default is `TRUE`.
#' @param method A character string specifying the correlation method (`"spearman"`, `"pearson"`, or `"kendall"`). Default is `"spearman"`.
#' @param r.threshold A numeric value specifying the minimum correlation coefficient for displaying results. Default is `0`.
#' @param p.threshold A numeric value specifying the maximum p-value for displaying results. Default is `0`.
#' @param theme.size An integer specifying the text size for axis labels and legends. Default is `10`.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item `p1`: A heatmap visualizing the correlation results.
#'   \item `p2`: A scatterplot-style visualization of the correlation results.
#'   \item `top10`: A matrix of the top OTUs/taxa retained for the analysis.
#'   \item `data`: A long-format data frame of the correlation results used for plotting.
#' }
#'
#' @examples
#' \dontrun{
#' result <- cor_env_ggcorplot_top(ps = ps, env1 = env1, Top = 5, method = "spearman")
#'
#' # Visualize the results
#' print(result$p1)  # Heatmap
#' print(result$p2)  # Scatterplot-style visualization
#' }
#'
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @export
cor_env_ggcorplot_top = function(ps =  ps.tem,
                                 jj = NULL,
                                 tran  =TRUE,
                                 Top = 10,
                                 env1 = ftab,
                                 label =  TRUE,
                                 col_cluster = TRUE,
                                 row_cluster = TRUE,
                                 method = "spearman",
                                 r.threshold= 0,
                                 p.threshold= 0,

                                 theme.size = 10){

  if (is.null(jj)) {
    psdata <- ps
  } else{
    psdata <- ggClusterNet::tax_glom_wt(ps = ps,ranks = jj)
  }

  if (tran) {
    psdata = psdata %>%
      phyloseq::transform_sample_counts(function(x) {x/sum(x)})
  }

  otu = phyloseq::otu_table(psdata)
  tax = phyloseq::tax_table(psdata)

  if (dim(otu)[1] < Top) {
    top10 <- otu[names(sort(rowSums(otu), decreasing = TRUE)[1:dim(otu)[1]]),]
    top10 = t(top10)
  } else {
    top10 <- otu[names(sort(rowSums(otu), decreasing = TRUE)[1:Top]),]
    top10 = t(top10)
  }
  head(top10)


  env2 = top10
  env2 <- env2[match(row.names(env1),row.names(env2)),]
  env0 <- cbind(env1,env2)
  occor = psych::corr.test(env0,use="pairwise",method=method,adjust="fdr",alpha=.05)
  occor.r = occor$r
  occor.p = occor$p

  occor.r[occor.p > p.threshold&abs(occor.r) < r.threshold] = 0

  head(env0)
  # data[data > 0.3]<-0.3
  #drop gene column as now in rows

  if (col_cluster ==  TRUE) {
    clust <- hclust(dist(env1 %>% as.matrix()%>% t())) # hclust with distance matrix
    ggtree_plot <- ggtree::ggtree(clust)
  }
  if (row_cluster ==  TRUE) {
    v_clust <- hclust(dist(env2 %>% as.matrix() %>% t()))
    ggtree_plot_col <- ggtree::ggtree(v_clust) + ggtree::layout_dendrogram()
  }

  head(occor.r)
  occor.r = as.data.frame(occor.r)
  data <- occor.r[colnames(env1),colnames(env2)]
  data$id = row.names(data)
  pcm = reshape2::melt(data, id = c("id"))

  head(pcm)

  p1 = ggplot(pcm, aes(y = id, x = variable)) +
    # geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) +
    geom_tile(aes(size = value,fill = value), guide = FALSE)+
    scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1), guide = FALSE) +
    labs( y= "", x = "",
          # size = "Relative Abundance (%)",
          fill = "")  +
    # scale_fill_manual(values = colours, guide = FALSE) +
    scale_x_discrete(limits = rev(levels(pcm$variable)))  +
    scale_y_discrete(position = "right") +
    scale_fill_gradientn(colours =colorRampPalette(c("#377EB8","#F7F4F9","#E41A1C"))(60)) +
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.x = element_text(colour = "black",size = theme.size,angle = 60,vjust = 1,hjust = 1)
    )

  # colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
  #              "#F6AE2D","#86BBD8")
  #----样本在y轴上
  p2 = ggplot(pcm, aes(y = id, x = variable)) +
    geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) +
    scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) +
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  +
    # scale_fill_manual(values = colours, guide = FALSE) +
    scale_x_discrete(limits = rev(levels(pcm$variable)))  +
    scale_y_discrete(position = "right")  +
    scale_fill_gradientn(colours =colorRampPalette(c("#377EB8","#F7F4F9","#E41A1C"))(60))  +
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.x = element_text(colour = "black",size = theme.size,angle = 60,vjust = 1,hjust = 1)
    )

  if (col_cluster ==  TRUE) {
    p1 <- p1  %>%
      aplot::insert_left(ggtree_plot, width=.2)
    p2 <- p2  %>%
      aplot::insert_left(ggtree_plot, width=.2)
  }

  if (label ==  TRUE) {
    p1 <- p1  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
    p2 <- p2  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
  }


  return(list(p1,p2,top10,data))


}
