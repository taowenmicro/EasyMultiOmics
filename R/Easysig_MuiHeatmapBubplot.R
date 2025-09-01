#' MuiHeatmapBubplot
#'
#' @param data a data.frame contain the input data
#' @param  i col index wtich need to test
#' @param errbar add a distinctive label, TRUE or FEASE could be selected
#' @param result output from aovMcomper or KwWlx. You can also import result calculated from other software (a data frame)
#' @param sample sample or group
#' @param col_cluster only used for sample = TRUE
#' @param row_cluster only used for sample = TRUE
#' @param label only used for sample = TRUE
#' @examples
#' # data(data_wt)
#'result = MuiKwWlx(data = data_wt,num = c(4,9,8))
#'res <- MuiHeatmapBubplot(data = data_wt,i = c(4,9,8) ,result = result,
#'                      sample  = TRUE)
#'res[[1]]
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} PengHao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


MuiHeatmapBubplot = function(
  data = data,
  i,
  result = result,
  sample = FALSE,
  scale = TRUE,
  col_cluster = TRUE,
  row_cluster = TRUE,
  label = TRUE
){
  mean = function(x) {
    base::mean(x,na.rm = TRUE)
  }
  i = c(1,2,i)
  data <- data[i]
  # data[data > 0.3]<-0.3
  mat <- data[,c(-1,-2)] #drop gene column as now in rows
  if (col_cluster ==  TRUE) {
    row.names(mat) = data$ID
    clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix
    ggtree_plot <- ggtree::ggtree(clust)
  }
  if (row_cluster ==  TRUE) {
    v_clust <- hclust(dist(mat %>% as.matrix() %>% t()))
    ggtree_plot_col <- ggtree::ggtree(v_clust) + ggtree::layout_dendrogram()
  }

  if (label ==  TRUE) {
    map = data[,c(1,2)]
    colnames(map)[1] = "ID"
    labels= ggplot(map, aes(x = 0, y=ID, fill=group)) + geom_tile() +
      scale_fill_brewer(palette = 'Set1',name="Cell Type") +
      theme_void()
  }

  if (sample == TRUE) {
    tem = data[,c(-2)]
    colnames(tem)[1] = "id"

    if (scale == TRUE) {
      tem.2 = scale((tem[,c(-1)])) %>%
        as.data.frame()
      tem = cbind(id = tem[,c(1)],tem.2)
    }

    pcm = reshape2::melt(tem, id = c("id"))
    head(pcm)
    # pcm$variable = factor(pcm$variable,levels = map$ID)
    # pcm$id = factor(pcm$id,levels = rig$id)

    p1 = ggplot(pcm, aes(y = id, x = variable)) +
      # geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) +
      geom_tile(aes(fill = value))+
      scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) +
      labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  +
      # scale_fill_manual(values = colours, guide = FALSE) +
      scale_x_discrete(limits = rev(levels(pcm$variable)))  +
      scale_y_discrete(position = "right") +
      scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral")[11:1])(60))+
      theme(
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.text.x = element_text(colour = "black",angle = 90,hjust = 1,
                                   vjust = 0)

      )

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
        axis.text.x = element_text(colour = "black",angle = 90,hjust = 1,
                                   vjust = 0)

      )

    if (label ==  T) {
      p1 <- p1  %>%
        aplot::insert_left(labels, width=.05)
      p2 <- p2  %>%
        aplot::insert_left(labels, width=.05)
    }


    if (col_cluster ==  T) {
      p1 <- p1  %>%
        aplot::insert_left(ggtree_plot, width=.2)
      p2 <- p2  %>%
        aplot::insert_left(ggtree_plot, width=.2)
    }



    if (row_cluster ==  T) {
      p1 <- p1  %>%
        aplot::insert_top(ggtree_plot_col, height=.2)
      p2 <- p2  %>%
        aplot::insert_top(ggtree_plot_col, height=.2)
    }

    dat = pcm
  } else {
    tem = data[,c(-1)]
    colnames(tem)[1] = "id"

    if (scale == TRUE) {
      tem.2 = scale((tem[,c(-1)])) %>%
        as.data.frame()
      tem = cbind(id = tem[,c(1)],tem.2)
    }

    result$group = row.names(result)
    abc <- reshape2::melt(result, id="group", variable.name="variable", value.name = "abc")
    colnames(abc)[1] = "id"
    tem.2 = tem %>% dplyr::group_by(id) %>%
      dplyr::summarise_if(is.numeric, mean) %>%
      as.data.frame()

    pcm = reshape2::melt(tem.2, id = c("id"))
    head(abc)

    dat = pcm %>% dplyr::left_join(abc)

    p1 = ggplot(dat, aes(x = id, y = variable)) +
      # geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) +
      geom_tile(aes(fill = value)) +
      geom_text(aes(label = abc)) +
      scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral")[11:1])(60))+
      theme(
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.text.x = element_text(colour = "black",angle = 90)
      )

    #----样本在y轴上
    p2 = ggplot(dat, aes(x = id, y = variable)) +
      geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) +
      geom_text(aes(label = abc)) +
      scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) +
      labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  +
      # scale_y_discrete(position = "right")  +
      scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral")[11:1])(60)) +
      theme(
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.text.x = element_text(colour = "black",angle = 90)

      )


  }
  return(list(p1,p2,dat))
}









