#' @title Generate Bar and Flow Plots for Microbial Composition
#'
#' @description
#' This function generates bar plots, alluvial flow plots, and flower plots to visualize microbial composition across groups.
#' The function supports relative abundance transformations, error bars, and optional labeling.
#'
#' @param otu A data frame containing OTU data. If `NULL`, the `ps` object is used.
#' @param tax A data frame containing taxonomy data. If `NULL`, the `ps` object is used.
#' @param map A data frame containing sample metadata. If `NULL`, the `ps` object is used.
#' @param tree A phylogenetic tree object. If `NULL`, the `ps` object is used.
#' @param ps A `phyloseq` object containing OTU, tax, mapping data.
#' @param group A character string specifying the grouping variable in the sample metadata. Default is `"Group"`.
#' @param j A character string specifying the taxonomic rank to analyze (e.g., `"Phylum"`, `"Genus"`). Default is `"Phylum"`.
#' @param axis_ord A character vector specifying the order of x-axis groups. If `NULL`, default order is used.
#' @param label Logical. If `TRUE`, adds text labels to the plots. Default is `TRUE`.
#' @param sd Logical. If `TRUE`, adds error bars to the plots. Default is `FALSE`.
#' @param Top An integer specifying the number of top taxa to display. Taxa not in the top are grouped as `"others"`. Default is `10`.
#' @param tran Logical. If `TRUE`, transforms data to relative abundances. Default is `TRUE`.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{Bar Plot}{A bar plot showing the relative abundances of taxa across groups.}
#'   \item{Processed Data}{A data frame with the processed microbial composition data.}
#'   \item{Alluvial Flow Plot}{An alluvial flow plot visualizing the changes in relative abundances across groups.}
#'   \item{Flower Plot}{A flower plot displaying relative abundances in a circular layout.}
#' }
#'
#' @details
#' The function processes microbiome data to compute relative abundances (if `tran = TRUE`) and groups taxa based on the specified rank (`j`).
#' It visualizes the data as bar plots, alluvial flow plots, and flower plots. The top `Top` taxa are highlighted, and the remaining taxa are grouped as `"others"`.
#'
#' @examples
#' \dontrun{
#' result = barMainplot.micro(ps = pst, j = "Genus", label = FALSE, sd = FALSE, Top = 5)
#' p4_1 <- result[[1]]
#' databar <- result[[2]]
#' }
#'
#' @author Contact: Tao Wen \email{2018203048@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njau.edu.cn}
#'
#' @export


# j = "Phylum"
# j = "Class"
# j = "Order"
# j =  "Family"
# j = "Genus"
#
# otu =NULL
# tax = NULL
# map = NULL
# ps = ps1
#
# group = "Group"
# axis_ord = NULL
# label = FALSE
# sd = FALSE
# Top = 10
#
# tran = TRUE

# library(phyloseq)
# library(tidyverse)
# library(vegan)
# library(reshape2)
# library("plyr")

# detach("package:ggalluvial")
# detach("package:plyr")
# detach("package:reshape2")
# library(ggplot2)


barMainplot.micro = function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = NULL,group  = "Group",
                       j = "Phylum",axis_ord = NULL,label = TRUE ,sd = FALSE,Top = 10,tran = TRUE){


  if (is.null(axis_ord)) {
    axis_order = NA
  }else{
    axis_order = strsplit(basename(axis_ord), "-")[[1]]
  }

  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)

  # psdata = ps %>%
  #   tax_glom(taxrank = j)
  psdata <- tax_glom_wt(ps = ps
                        ,ranks = j)

  # transform to relative abundance
  if (tran == TRUE) {

    psdata = psdata%>%
      phyloseq::transform_sample_counts(function(x) {x/sum(x,na.rm = TRUE)} )

  }


  otu = phyloseq::otu_table(psdata)
  tax = phyloseq::tax_table(psdata)
  # otu [1:20,1:20]
  dim(otu)
  otu[is.na(otu)] = 0
  for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {

      tax[i,j] =tax[i,j]
    } else {
      tax[i,j]= "others"
    }
  }
  phyloseq::tax_table(psdata)= tax

  Taxonomies <- psdata %>% # Transform to rel. abundance
    phyloseq::psmelt()


  if (tran == TRUE) {
    Taxonomies$Abundance = Taxonomies$Abundance * 100
  }
  # Taxonomies$Abundance = Taxonomies$Abundance /rep

  colnames(Taxonomies)[match(j,colnames(Taxonomies))] <- "aa"
  data = c()
  i = 2
  for (i in 1:length(unique(phyloseq::sample_data(ps)$Group))) {
    a <- as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,1]

    b =  as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,2]
    head(Taxonomies)
    colnames(Taxonomies)
    c <- Taxonomies %>%
      dplyr::filter(Group == as.character(a))
    c$Abundance <- c$Abundance/b
    data = data.frame(Sample =c$Sample,Abundance = c$Abundance,aa =c$aa,Group = c$Group)

    if (i == 1) {
      table = data
    }
    if (i != 1) {
      table = rbind(table,data)
    }

  }
  # head(table)
  # sum(table$Abundance)
  # Taxonomies$Abundance = NULL
  # Taxonomies <- Taxonomies %>% inner_join(table,by = "Sample")

  Taxonomies = table
  head(Taxonomies)

  #按照分组求均值

  by_cyl <- dplyr::group_by(Taxonomies, aa,Group)

  zhnagxu2 = dplyr::summarise(by_cyl, sum(Abundance,na.rm = TRUE), stats::sd(Abundance,na.rm = TRUE))
  head(zhnagxu2)

  iris_groups<- dplyr::group_by(Taxonomies, aa)
  cc<- dplyr::summarise(iris_groups, sum(Abundance,na.rm = TRUE))
  head(cc)
  colnames(cc)= c("aa","allsum")
  cc<- dplyr::arrange(cc, desc(allsum))

  ##使用这个属的因子对下面数据进行排序
  head(zhnagxu2)
  colnames(zhnagxu2) <- c("aa","group","Abundance","sd")
  zhnagxu2$aa = factor(zhnagxu2$aa,order = TRUE,levels = cc$aa)

  # zhnagxu3 = zhnagxu2
  zhnagxu3 = zhnagxu2 %>% left_join(cc) %>% arrange(allsum)


  # zhnagxu3 = zhnagxu3 %>% arrange(desc(aa))
  head(zhnagxu3)
  ##制作标签坐标，标签位于顶端
  # Taxonomies_x = ddply(zhnagxu3,"group", summarize, label_y = cumsum(Abundance))
  # head(Taxonomies_x )
  #标签位于中部
  # Taxonomies_x1 = ddply(zhnagxu3,"group", transform, label_y = cumsum(Abundance) - 0.5*Abundance)
  dat = plyr::ddply(zhnagxu3,'group',
                             dplyr::summarize,label_sd = cumsum(Abundance),
                             label_y = cumsum(Abundance) - 0.5*Abundance)


 # dat = zhnagxu3 %>% group_by(aa,group) %>%
 #    summarise(label_sd = cumsum(Abundance),
 #              label_y = cumsum(Abundance) - 0.5*Abundance
 #              )


  dat$aa = zhnagxu3$aa %>% unique()
  head(dat)
  # Taxonomies_x$group = NULL
  # Taxonomies_x = Taxonomies_x %>% arrange(desc(aa))
  # Taxonomies_x$aa = NULL
  # Taxonomies_x$label_y =

  Taxonomies_x = zhnagxu3 %>% inner_join(dat)


  Taxonomies_x$label = Taxonomies_x$aa
  # #使用循环将堆叠柱状图柱子比较窄的别写标签，仅仅宽柱子写上标签
  # for(i in 1:nrow(Taxonomies_x)){
  #   if(Taxonomies_x[i,3] > 3){
  #     Taxonomies_x[i,5] = Taxonomies_x[i,5]
  #   }else{
  #     Taxonomies_x[i,5] = NA
  #   }
  # }



  Taxonomies_x$aa = factor(Taxonomies_x$aa,order = TRUE,levels = c(as.character(cc$aa)))



  ##普通柱状图


  p4 <- ggplot(Taxonomies_x , aes(x =  group, y = Abundance, fill = aa, order = aa)) +
    geom_bar(stat = "identity",width = 0.5,color = "black") +
    # scale_fill_manual(values = colors) +
    theme(axis.title.x = element_blank()) +
    theme(legend.text=element_text(size=6)) +
    scale_y_continuous(name = "Relative abundance (%)") +
    guides(fill = guide_legend(title = j)) +
    labs(x="",y="Relative abundance (%)",
         title= "") +
    scale_fill_brewer(palette = "Paired")
  # paste("Top ",Top," ",j,sep = "")
  p4
  if (is.na(axis_order)) {
    p4 = p4
  }else{
    p4 = p4 + scale_x_discrete(limits = axis_order)
  }


  if (sd == TRUE) {
    p4 =  p4 +
      geom_errorbar(aes(ymin=label_sd-sd, ymax=label_sd +sd), width=.2)
  }

  if (label == TRUE) {
    p4 = p4 +
      geom_text(aes(y = label_y, label = label ))
  }


  map = as.data.frame(phyloseq::sample_data(ps))
  if (length(unique(map$Group))>3){p4 = p4 + theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))}

  #-------------冲击图
  cs = Taxonomies_x $aa

  lengthfactor <- cs %>%
    levels() %>%
    length()
  cs4 <- cs %>%
    as.factor() %>%
    summary()  %>%
    as.data.frame()
  cs4$id = row.names(cs4)


  #对因子进行排序
  df_arrange<- dplyr::arrange(cs4, id)
  #对Taxonomies_x 对应的列进行排序
  Taxonomies_x1<- dplyr::arrange(Taxonomies_x , aa)
  head(Taxonomies_x1)
  #构建flow的映射列Taxonomies_x
  Taxonomies_x1$ID = factor(rep(c(1:lengthfactor), cs4$.))

  #colour = "black",size = 2,,aes(color = "black",size = 0.8)
  head(Taxonomies_x1)
  Taxonomies_x1$Abundance

  p3 <- ggplot(Taxonomies_x1, aes(x = group, y = Abundance,fill = aa,alluvium = aa,stratum = ID)) +
    ggalluvial::geom_flow(aes(fill = aa, colour = aa),
              stat = "alluvium", lode.guidance = "rightleft",
              color = "black",size = 0.2,width = 0.35,alpha = .2)  +
    geom_bar(width = 0.45,stat = "identity") +
    labs(x="",y="Relative abundance (%)",
         title= "") +
    guides(fill = guide_legend(title = j),color = FALSE) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_brewer(palette = "Paired")
  p3


  # flower plot
  p1 <- ggplot(Taxonomies_x1,
               aes(x = group, alluvium = aa, y = Abundance)) +
    ggalluvial::geom_flow(aes(fill = aa, colour = aa), width = 0)  +
    labs(x="",y="Relative abundance (%)",
         title="") +
    guides(fill = guide_legend(title = j),color = FALSE) +
    scale_y_continuous(expand = c(0,0))


  if (is.na(axis_order)) {
    p1 = p1
    p3 = p3

  }else{
    p1 = p1 + scale_x_discrete(limits = axis_order)
    p3 = p3 + scale_x_discrete(limits = axis_order)

  }
  # p3
  if (label == TRUE) {
    p1 = p1 +
      geom_text(aes(y = label_y, label = label ))
    p3 = p3 +
      geom_text(aes(y = label_y, label = label ))
  }

  if (sd == TRUE) {
    p1 =  p1 +
      geom_errorbar(aes(ymin=label_sd-sd, ymax=label_sd +sd), width=.2)
    p3 =  p3 +
      geom_errorbar(aes(ymin=label_sd-sd, ymax=label_sd +sd), width=.2)
  }

  if (length(unique(map$Group))>3){	p3=p3+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))}

  return(list(p4,Taxonomies,p3,p1))
}




