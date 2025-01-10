#' @title  Create a heatmap with metabolites clustering and labeling
#'
#' @description
#' The heatmap.ms function generates a heatmap with clustering and labeling based on metabolites  data and sample metadata.
#'  It includes row and column clustering, relative abundance visualization, and sample group labeling.
#' @param ps_rela A phyloseq format file used as an alternative for the input containing metabolite composition table(realtive abundance),
#' metabolite classification table, and sample metadata.
#' @param label Logical parameters indicating whether to label sample groups. Default is TRUE.
#' @param col_cluster Logical parameters indicating whether to perform column clustering. Default is TRUE.
#' @param row_cluster Logical parameters indicating whether to perform row clustering. Default is TRUE.
#' @return A list containing the following components:
#' \item{p1}{Heatmap with row clustering and relative abundance visualization.}
#' \item{p2}{Bubble heatmap with sample clustering and relative abundance visualization.}
#' \item{pcm}{Analysis data for heatmap plotting.}
#'
#' @examples
#' ps.ms_rela <- ps.ms %>% scale_micro(method = "rela") %>%
#'  tax_glom_wt(ranks = "Class")
#'  result <- heatmap.ms (ps.ms_rela,label =  FALSE,col_cluster = FALSE,row_cluster =FALSE)
#'  p19 <- result[[1]]
#' @export

# ps_rela <- ps %>% scale_micro(method = "rela") %>%
#   tax_glom_wt(ranks = "classification")
#
# result <- GCheatmap (ps_rela,
#                       label =  F,
#                       col_cluster = F,
#                       row_cluster = F)
# p1 <- result[[1]]
# p1
# # p1 +  scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(60))
# p2 <- result[[2]]
# p2


heatmap.ms <- function(ps_rela= ps_rela,
                      label =  TRUE,
                      col_cluster = TRUE,
                      row_cluster = TRUE
){



  # print(ps_heatm)
  datah <- as.data.frame(t(vegan_otu(ps_rela)))
  head(datah)
  tax = as.data.frame(vegan_tax(ps_rela))
  otutaxh = cbind(datah,tax)
  otutaxh$id = row.names(otutaxh)
  # otutaxh$id  = row.names(otutaxh)
  row.names(otutaxh) = otutaxh$id
  data <- otutaxh[,c("id",sample_names(ps))]

  rig <- data[,sample_names(ps_rela)] %>% rowMeans() %>% as.data.frame()
  head(rig)
  colnames(rig) = "MeanAbundance"
  rig$id = row.names(rig)
  rig = rig %>% arrange(MeanAbundance)
  rig$id = factor(rig$id,levels = rig$id)
  p_rig = ggplot(rig) + geom_bar(aes(y = id,x = MeanAbundance),
                                 fill = "#A54657",
                                 stat = "identity") + theme_void()

  tem = data[,sample_names(ps)] %>% as.matrix()

  tem = scale(t(tem)) %>% t() %>%
    as.data.frame()
  data[,sample_names(ps)] = tem



  # data[data > 0.3]<-0.3
  mat <- data[,-1] #drop gene column as now in rows

  if (col_cluster ==  TRUE) {
    clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix
    ggtree_plot <- ggtree::ggtree(clust)
  }

  if (row_cluster ==  TRUE) {
    v_clust <- hclust(dist(mat %>% as.matrix() %>% t()))
    ggtree_plot_col <- ggtree::ggtree(v_clust) + ggtree::layout_dendrogram()
  }

  if (label ==  TRUE) {
    map = as.data.frame(sample_data(ps))
    map$ID = row.names(map)
    labels= ggplot(map, aes(x = ID, y=1, fill=Group)) + geom_tile() +
      scale_fill_brewer(palette = 'Set1',name="Cell Type") +
      theme_void()
  }

  map = sample_data(ps_rela) %>% as.tibble() %>%
    arrange(Group) %>% as.data.frame()
  map$ID
  pcm = reshape2::melt(data, id = c("id"))
  pcm$variable = factor(pcm$variable,levels = map$ID)
  pcm$id = factor(pcm$id,levels = rig$id)


  p1 = ggplot(pcm, aes(y = id, x = variable)) +
    # geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) +
    geom_tile(aes(size = value,fill = value))+
    scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) +
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  +
    # scale_fill_manual(values = colours, guide = FALSE) +
    scale_x_discrete(limits = rev(levels(pcm$variable)))  +
    scale_y_discrete(position = "right") +
    scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral")[11:1])(60))+
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.x = element_text(colour = "black",angle = 90)

    )

  colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
               "#F6AE2D","#86BBD8")
  #----样本在y轴上
  p2 = ggplot(pcm, aes(y = id, x = variable)) +
    geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) +
    scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) +
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  +
    # scale_fill_manual(values = colours, guide = FALSE) +
    scale_x_discrete(limits = rev(levels(pcm$variable)))  +
    scale_y_discrete(position = "right")  +
    scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral")[11:1])(60)) +
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.x = element_text(colour = "black",angle = 90)

    )

  p1 <- p1  %>%
    aplot::insert_right(p_rig, width=.2)
  p2 <- p2  %>%
    aplot::insert_right(p_rig, width=.2)
  if (col_cluster ==  TRUE) {
    p1 <- p1  %>%
      aplot::insert_left(ggtree_plot, width=.2)
    p2 <- p2  %>%
      aplot::insert_left(ggtree_plot, width=.2)
  }

  if (row_cluster ==  TRUE) {
    p1 <- p1  %>%
      aplot::insert_top(labels, height=.02)
    p2 <- p2  %>%
      aplot::insert_top(labels, height=.02)
  }

  if (label ==  TRUE) {
    p1 <- p1  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
    p2 <- p2  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
  }

  return(list(p1,p2,pcm))
}


