
#' @title Paired Volcano Plot for Correlation Analysis Across Different Datasets
#'
#' @param ps01 A `phyloseq` object containing microbiome data.
#' @param ps02 A `phyloseq` object containing metabolomics data.
#' @param group1 Character. Label for the first dataset (e.g., `"micro"` for microbiome).
#' @param group2 Character. Label for the second dataset (e.g., `"ms"` for metabolomics).
#' @param r.threshold Numeric. Correlation coefficient threshold to include in analysis (default is 0.6).
#' @param p.threshold Numeric. P-value threshold for significance (default is 0.1).
#' @param method Character. Correlation method (`"spearman"`, `"pearson"`, etc.).
#' @param top Integer. Number of top features to consider for analysis (default is 500).
#'
#' @return A list containing:
#' \itemize{
#'   \item Plots: Volcano plots with significant features and correlation connections.
#'   \item Data: Dataframes of significant features and their statistics.
#'   \item Edges: Dataframes of correlations between features.
#' }
#' @examples
#' \dontrun{
#'  results <- volcano.line2.omics(ps01,ps02, group1 = "micro", group2 = "ms")
#'
#' }
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}

volcano.line2.omics<- function(ps01=ps01,
                               ps02 =ps02,
                               group1= "micro",
                               group2= "ms",
                               r.threshold = 0.6,
                               p.threshold = 0.1,
                               method = "spearman",
                               top = 500){

  # 计算差异微生物

  res = EdgerSuper.micro(ps = ps01 %>% ggClusterNet::filter_OTU_ps(top),group  = "Group",
                         artGroup = NULL, j = "OTU")

  da16 = res[[2]]
  da16$id =  row.names(da16)
  row.names(da16)= NULL
  # 差异代谢物
  dams = statSuper(ps =ps02,group  = "Group",artGroup = NULL,method = "wilcox")
  head(dams)
  head(da16)
  # 提取分组
  sub_design <- as.data.frame(phyloseq::sample_data(ps01))
  Desep_group <- as.character(levels(as.factor(sub_design$Group)))

  aaa = combn(Desep_group,2)

  # results <- list()
  plot1 = list()
  data1 = list()
  data2 = list()

  for (ii in 1:dim(aaa)[2]) {
    # ii = 1
    Desep_group = aaa[,ii]
    print( Desep_group)

    col_16 <- paste(paste(Desep_group[1], Desep_group[2], sep = "-"), c("logFC", "p"), sep = "")

    # 筛选ko富集的
    da16_1 =  da16 %>% dplyr::select(col_16,id) %>%
      dplyr::filter( .data[[col_16[1]]] >0 ) %>% mutate(group1=group1)
    head(da16_1)
    da16_1[is.na(da16_1)] = 0


    col_ms <- paste(Desep_group[1], Desep_group[2], c("log2_FC", "Pvalue"), sep = "_")
    head(dams)
    # 筛选ko富集的
    dams_1 =  dams%>% dplyr::select(col_ms,id) %>% dplyr::filter( .data[[col_ms[1]]] >0 ) %>% mutate(group1=group2)
    dams_1[is.na(dams_1)] = 0
    dams_1[[col_ms[1]]] [dams_1[[col_ms[1]]] == "NA"] = 0
    dams_1[[col_ms[1]]] = as.numeric(dams_1[[col_ms[1]]])
    dams_1[[col_ms[1]]]=  - (dams_1[[col_ms[1]]])
    head(dams_1)

    dams_1[[col_ms[2]]] [dams_1[[col_ms[2]]] == "NA"] = 1
    dams_1[[col_ms[2]]] = as.numeric(dams_1[[col_ms[2]]])
    colnames(dams_1) =   colnames(da16_1)

    dat = rbind(dams_1,da16_1)
    colnames(dat)[1:2]= c("logFC","pvalue")

    dat <- dat %>%
      mutate(level = ifelse(pvalue < 0.05, "enrich", "nosig"))
    dat $level2 = paste(dat $level, dat $group1 ,sep = "_")

    dat <- dat %>%
      mutate(level2 = ifelse(str_detect(level2,"nosig") , "nosig",level2))

    head(dat)

    p= ggplot( dat,aes(x =logFC ,y = -log2(pvalue)
                       ,fill=level2
    )) +
      geom_point(aes(shape = level2 ),size =2,color = "black" )+
      scale_shape_manual(values = c( 22,24,1)) +
      theme_bw()+
      scale_fill_manual(values = c("#ffad73","#26b3ff","grey90"))

    p

    # 添加连线

    # 分别提取对应的ps文件
    sample_data(ps01)

    ps01_2 =  ps01 %>%  subset_samples.wt("Group",Desep_group ) %>%
      subset_taxa.wt("OTU", row.names( phyloseq::tax_table(ps01)) %in% dat$id)

    ps02_2 = ps02 %>%   subset_samples.wt("Group",Desep_group )  %>%
      subset_taxa.wt("OTU", row.names( phyloseq::tax_table(ps02)) %in% dat$id)


    ps_sub = merge.ps(ps01_2,ps02_2,onlygroup = TRUE,dat1.lab = group1,  dat2.lab = group2 )

    result <- corBiostripeBig(ps = ps_sub, r.threshold = r.threshold,
                              p.threshold = p.threshold, method = method)

    cor = result [[1]]

    colnames(cor) = gsub(paste(group1,"_",sep = ""),"",colnames(cor))
    colnames(cor) = gsub(paste(group2,"_",sep = ""),"",colnames(cor))


    row.names(cor) = gsub(paste(group1,"_",sep = ""),"",row.names(cor))
    row.names(cor) = gsub(paste(group2,"_",sep = ""),"",row.names(cor))


    #
    # id1= colnames(cor) %in% (dat %>% dplyr::filter(group1==group2) %>% .$id)
    # id2= colnames(cor) %in% (dat %>% dplyr::filter(group1!=group2) %>% .$id)
    #
    # cor2 = cor[id2,id1]

    res = dat

    node = data.frame(elements =res$id,X1 = res$logFC,
                      X2 =res$pvalue)

    edge = edgeBuild(cor = cor,node = node)

    head(edge)


    # ggplot() +
    #   geom_point(dat,aes(x = logFC , y= -log2(pvalue) ,fill=level2,shape = level2 ),
    #              size =2, color = "black" ) +
    #   scale_shape_manual(values = c( 22,24,1)) +
    #   theme_bw()+
    #   scale_fill_manual(values = c("#ffad73","#26b3ff","grey90"))+
    #   ggnewscale::new_scale_color() +
    #   geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
    #                data = edge, size = 0.5)+
    #   scale_color_manual(values = c('blue2','red2', 'gray30'))+ggtitle(paste(group1, group2,sep=""))
    #
    #
    # ggplot(dat, aes(x = logFC, y = -log2(pvalue), fill = level2, shape = level2)) +
    #   geom_point(size = 2, color = "black") +
    #   scale_shape_manual(values = c(22, 24, 1)) +
    #   scale_fill_manual(values = c("#ffad73", "#26b3ff", "grey90")) +
    #   theme_bw() +
    #   ggnewscale::new_scale_color() +
    #   geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2, color = as.factor(cor)),
    #                data = edge, size = 0.5) +
    #    scale_color_manual(values = c('blue2', 'red2', 'gray30')) +
    #   ggtitle(paste(group1, group2, sep = ""))


    dat2 = dat %>% dplyr::filter(level!= "nosig")

    edge2 =  edge %>% dplyr::filter( OTU_2%in% dat2 $id |OTU_1%in% dat2 $id   )

    dat3 = dat2 %>% dplyr::filter( id%in% c(edge2 $OTU_2,edge2 $OTU_1)  )


    p1= ggplot()+geom_segment(aes(x = X1, y = -log2(Y1), xend = X2, yend = -log2(Y2), color = as.factor(cor)),
                              data = edge2, size = 0.2,alpha = 0.5, linetype = "dashed") +
      scale_color_manual(values = c('blue2', 'red2')) +
      ggnewscale::new_scale_color() +
      geom_point(data = dat,aes(x = logFC , y= -log2(pvalue) , shape = level2,
                                #      group = level2
                                fill=level2  ), alpha = 0.7)+
      ggrepel::geom_text_repel(data =dat3,aes(x = logFC , y= -log2(pvalue), label = id))+
      scale_fill_manual(values = c("#ffad73", "#26b3ff", "grey90")) +
      scale_shape_manual(values = c(22, 24, 1))+
      ggtitle(paste(Desep_group, collapse = "_"))+
      theme_bw()+
      labs(   x = "Log2 fold change",    # 修改横坐标标签
              y = "-Log10 P",
              color = "Cor")
    p1

    print("1")

    plot1[[ii]]=  p1

    names(plot1)[ii] <- paste(Desep_group, collapse = "_")


    data1[[ii]]=  dat

    names(data1)[ii] <- paste(Desep_group, collapse = "_")

    data2[[ii]]=  edge
    names(data2)[ii] <- paste(Desep_group, collapse = "_")


    #
    # results[[paste(Desep_group, collapse = "_")]] <- dat
    #
    # results[[paste(Desep_group, collapse = "_")]] <- p1
    #
    # results[[paste(Desep_group, collapse = "_")]] <- edge

  }


  return(list(plot1,data1,data2))

}


