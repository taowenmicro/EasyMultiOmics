
library(ggalluvial)

tax_glom_trans <- function(ps = ps,ranks = "pathway") {


  if (  is.numeric(ranks)) {
    ranks <- phyloseq::rank.names(ps)[ranks]
  }


  otu <- as.data.frame(t(vegan_otu(ps)))
  tax <- as.data.frame(vegan_tax(ps))
  head(tax)
  # building group
  tax[[ranks]][is.na(tax[[ranks]])] = "Unknown"
  tax[[ranks]][tax[[ranks]] == ""] = "Unknown"
  tax[[ranks]][tax[[ranks]] == "NA"] = "Unknown"
  split <- split(otu,tax[[ranks]])
  #calculate sum by group
  apply <- lapply(split,function(x)colSums(x[]))
  # result chack
  otucon <- do.call(rbind,apply)

  head(tax)
  # row.names(otucon) %>% unique()
  # tax$level2 %>% unique()
  # tax[[ranks]] %>% unique()

  taxcon <- tax[1:match(ranks,colnames(tax))] %>% as.matrix()
  taxcon <- taxcon[!duplicated(tax[[ranks]]),] %>% as.data.frame()
  head(taxcon)
  # colnames(taxcon) = j
  if (ncol(taxcon) == 1) {
    colnames(taxcon) = j
  }
  #-tax name with NA wound be repeated with unknown
  taxcon[[ranks]][is.na(taxcon[[ranks]])] = "unknown"
  row.names(taxcon) <- taxcon[[ranks]]

  match(row.names(otucon),row.names(taxcon))
  head(tax)
  pscon <- phyloseq::phyloseq(
    phyloseq::otu_table( as.matrix(otucon),taxa_are_rows = TRUE),
    phyloseq::tax_table(as.matrix(taxcon)),
    phyloseq::sample_data(ps)
  )
  return(pscon)
}
#' @title Generate Bar and Flow Plots for Transcriptome Functional Composition Data
#'
#' @description
#' This function generates bar plots, alluvial flow plots, and flower plots to visualize transcriptome functional composition data across groups.
#' The function supports relative abundance transformations, error bars, and optional labeling.
#'
#' @param otu A data frame containing transcriptome functional composition table. If `NULL`, the `ps` object is used.
#' @param tax A data frame containing taxonomy data. If `NULL`, the `ps` object is used.
#' @param map A data frame containing sample metadata. If `NULL`, the `ps` object is used.
#' @param tree A phylogenetic tree object. If `NULL`, the `ps` object is used.
#' @param ps A phyloseq format file used as an alternative for the input containing transcriptome functional composition table, tax, and sample metadata.
#' @param group A character string specifying the grouping variable in the sample metadata. Default is `"Group"`.
#' @param j A character string specifying the taxonomic rank to analyze (e.g., `"Level1"`, `"Pathway"`). Default is `"Pathway"`.
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
#'   \item{Processed Data}{A data frame with the processed transcriptome functional composition data.}
#'   \item{Alluvial Flow Plot}{An alluvial flow plot visualizing the changes in relative abundances across groups.}
#'   \item{Flower Plot}{A flower plot displaying relative abundances in a circular layout.}
#' }
#'
#' @details
#' The function processes transcriptome functional composition data to compute relative abundances (if `tran = TRUE`) and groups taxa based on the specified rank (`j`).
#' It visualizes the data as bar plots, alluvial flow plots, and flower plots. The top `Top` taxa are highlighted, and the remaining taxa are grouped as `"others"`.
#'
#' @examples
#' \dontrun{
#' library(phyloseq)
#' library(dplyr)
#' library(ggplot2)
#' library(tidyverse)
#' library(ggClusterNet)
#' library(RColorBrewer)
#' # Example with a phyloseq object
#' data(ps.trans)
#' ps = ps.trans %>% filter_taxa(function(x) sum(x) > 5, TRUE)
#' result = barMainplot.trans(ps = ps, j = "level2", label = FALSE, sd = FALSE, Top = 10)
#' p4_1 <- result[[1]]
#' p4_2 <- result[[3]] + scale_fill_brewer(palette = "Paired")
#' databar <- result[[2]]
#' head(databar)
#' }
#'
#' @author Contact: Tao Wen \email{2018203048@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njau.edu.cn}
#'
#' @export


barMainplot.trans = function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = NULL,group  = "Group",
                            j = "Phylum",axis_ord = NULL,label = TRUE ,sd = FALSE,Top = 10,tran = TRUE){


  if (is.null(axis_ord)) {
    axis_order = NA
  }else{
    axis_order = strsplit(basename(axis_ord), "-")[[1]]
  }

  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)

  # psdata = ps %>%
  #   tax_glom(taxrank = j)

  psdata <- tax_glom_trans(ps = ps,ranks = j)

  # transform to relative abundance
  if (tran == TRUE) {
    psdata = psdata%>%
      phyloseq::transform_sample_counts(function(x) {x/sum(x)} )
  }


  otu = phyloseq::otu_table(psdata)
  tax = phyloseq::tax_table(psdata)

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


  Taxonomies$Abundance = Taxonomies$Abundance * 100
  # Taxonomies$Abundance = Taxonomies$Abundance /rep
  colnames(Taxonomies) <- gsub(j,"aa",colnames(Taxonomies))
  data = c()

  for (i in 1:length(unique(phyloseq::sample_data(ps)$Group))) {
    a <- as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,1]

    b =  as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,2]

    c <- Taxonomies %>%
      dplyr::filter(Group == a)
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
  # head(Taxonomies)
  Taxonomies = table



  #按照分组求均值

  by_cyl <- dplyr::group_by(Taxonomies, aa,Group)

  zhnagxu2 = dplyr::summarise(by_cyl, sum(Abundance), sd(Abundance))


  iris_groups<- dplyr::group_by(Taxonomies, aa)
  cc<- dplyr::summarise(iris_groups, sum(Abundance))
  head(cc)
  colnames(cc)= c("aa","allsum")
  cc<- dplyr::arrange(cc, desc(allsum))

  ##使用这个属的因子对下面数据进行排序
  head(zhnagxu2)
  colnames(zhnagxu2) <- c("aa","group","Abundance","sd")
  zhnagxu2$aa = factor(zhnagxu2$aa,order = TRUE,levels = cc$aa)

  zhnagxu3 = zhnagxu2

  ##制作标签坐标，标签位于顶端
  # Taxonomies_x = ddply(zhnagxu3,"group", summarize, label_y = cumsum(Abundance))
  # head(Taxonomies_x )
  #标签位于中部
  # Taxonomies_x1 = ddply(zhnagxu3,"group", transform, label_y = cumsum(Abundance) - 0.5*Abundance)
  Taxonomies_x = plyr::ddply(zhnagxu3,"group", dplyr::summarize,label_sd = cumsum(Abundance),label_y = cumsum(Abundance) - 0.5*Abundance)
  head( Taxonomies_x )
  zhnagxu3$group = as.factor(zhnagxu3$group)
  # Taxonomies_x$label_y =
  Taxonomies_x = cbind(as.data.frame(zhnagxu3),as.data.frame(Taxonomies_x)[,-1])


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
         title= "")
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
    scale_y_continuous(expand = c(0,0))
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



#
# tax_glom_wt <- function(ps = ps,ranks = "Phylum") {
#
#
#   otu <- as.data.frame(t(vegan_otu(ps)))
#   tax <- as.data.frame(vegan_tax(ps))
#
#   # building group
#   tax[[ranks]][is.na(tax[[ranks]])] = "Unknown"
#   tax[[ranks]][tax[[ranks]] == ""] = "Unknown"
#   split <- split(otu,tax[[ranks]])
#   #calculate sum by group
#   apply <- lapply(split,function(x)colSums(x[]))
#   # result chack
#   otucon <- do.call(rbind,apply)
#
#   taxcon <- tax[1:match(ranks,colnames(tax))]
#   taxcon <- taxcon[!duplicated(tax[[ranks]]),]
#   #-tax name with NA wound be repeated with unknown
#   taxcon[[ranks]][is.na(taxcon[[ranks]])] = "unknown"
#   row.names(taxcon) <- taxcon[[ranks]]
#
#
#   pscon <- phyloseq(
#     otu_table( as.matrix(otucon),taxa_are_rows = TRUE),
#     tax_table(as.matrix(taxcon)),
#     sample_data(ps)
#   )
#   return(pscon)
# }
#
