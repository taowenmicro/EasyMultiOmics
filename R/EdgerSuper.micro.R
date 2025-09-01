#' @title Differential Abundance Analysis Using EdgeR for Microbial Data
#'
#' @description
#' This function performs differential abundance analysis for microbial data using the
#' \code{EdgeR} package. It identifies differentially abundant taxa between groups in
#' a \code{phyloseq} object, visualizes the results with volcano plots, and outputs
#' normalized abundance data with taxonomic annotations.
#'
#' @param otu A data frame containing OTU counts. Optional if \code{ps} is provided.
#' @param tax A data frame containing taxonomic annotations. Optional if \code{ps} is provided.
#' @param map A data frame containing sample metadata. Optional if \code{ps} is provided.
#' @param tree A phylogenetic tree. Optional if \code{ps} is provided.
#' @param ps A \code{phyloseq} object containing microbiome data. If provided,
#'   overrides \code{otu}, \code{tax}, \code{map}, and \code{tree}.
#' @param group A string specifying the grouping variable in the sample metadata.
#'   Default is \code{"Group"}.
#' @param pvalue A numeric value specifying the significance threshold for adjusted
#'   p-values. Default is \code{0.05}.
#' @param lfc A numeric value specifying the log2 fold change cutoff. Default is \code{0}.
#' @param artGroup A matrix specifying custom group comparisons. Optional.
#' @param method A string specifying the normalization method for EdgeR. Default is \code{"TMM"}.
#' @param j A string or integer specifying the taxonomic rank to perform the analysis on.
#'   Can be a numeric rank (1-7), taxonomic name (e.g., \code{"Phylum"}), or \code{"OTU"}.
#'   Default is \code{2}.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{Plots}{A list of volcano plots for each group comparison, showing log fold
#'     changes vs. p-values.}
#'   \item{Results}{A data frame containing differential abundance results, including
#'     normalized abundances and taxonomic annotations.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Prepares OTU data, taxonomic annotations, and metadata from a \code{phyloseq}
#'     object or individual input files.
#'   \item Aggregates data to the specified taxonomic rank (if applicable).
#'   \item Normalizes OTU counts using the specified normalization method (\code{"TMM"}
#'     by default).
#'   \item Constructs a model matrix for group comparisons and estimates dispersions
#'     using EdgeR's generalized linear model (GLM) framework.
#'   \item Identifies differentially abundant taxa based on the specified \code{pvalue}
#'     and \code{lfc} thresholds.
#'   \item Generates volcano plots for each group comparison, highlighting enriched
#'     and depleted taxa.
#'   \item Appends taxonomic annotations and normalized abundances to the differential
#'     abundance results.
#' }
#'
#' The function supports both predefined and custom group comparisons (\code{artGroup})
#' and provides detailed visualizations for the results.
#'
#' @examples
#' \dontrun{
#' res <- EdgerSuper.micro(
#'   ps = ps.16s %>% ggClusterNet::filter_OTU_ps(500),
#'   group = "Group",
#'   artGroup = NULL,
#'   j = "OTU"
#' )
#' p25.1 <- res[[1]][1]
#' p25.2 <- res[[1]][2]
#' p25.3 <- res[[1]][3]
#' dat <- res[[2]]
#' }
#'
#' @author Contact: Tao Wen \email{2018203048@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njau.edu.cn}
#' @export




#------许多情况下，我们需要将指定一些组合，并不是全部组别的两两组合，然后做差异分析#--------
#-----输出我希望除了指定文件夹之外，可以将指定的分组全部组合到一起，然后输出到一张表格#---
# #--------差异分析:Edger#---------
#
# #--根据分组，将分组内内全部的组合都会做差异分析，并输出每个两两比较的csv表格（全部otu和差异otu表格）


# otu = NULL
# tax = NULL
# map = NULL
# tree = NULL
# ps = ps
# group  = "Group"
# pvalue = 0.05
# lfc =0
# artGroup = NULL

# method could be selected:  TMM,RLE, upperquartile.

EdgerSuper.micro = function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,
                      ps = NULL,group  = "Group",pvalue = 0.05,
                      lfc =0,artGroup = NULL,
                      method = "TMM",
                      j = 2
                      ){

  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)
  ps
  # ps = ps %>%
  #   ggClusterNet::tax_glom_wt(ranks = j)
  if (j %in% c("OTU","gene","meta")) {
    ps = ps
  } else if (j %in% c(1:7)) {
    ps = ps %>%
      ggClusterNet::tax_glom_wt(ranks = j)
  } else if (j %in% c("Kingdom","Phylum","Class","Order","Family","Genus","Species")){
    ps = ps %>%
      ggClusterNet::tax_glom_wt(ranks = j)
  } else {
    ps = ps
    print("unknown j, checked please")
  }
  sub_design <- as.data.frame(phyloseq::sample_data(ps))
  colnames(sub_design) = "Group"
  Desep_group <- as.character(levels(as.factor(sub_design$Group)))
  Desep_group


  if ( is.null(artGroup)) {
    #--构造两两组合#-----
    aaa = combn(Desep_group,2)
    # sub_design <- as.data.frame(sample_data(ps))
  }
  if (!is.null(artGroup)) {
    aaa  = as.matrix(b )
  }
  otu_table = as.data.frame(ggClusterNet::vegan_otu(ps))
  count = as.matrix(otu_table)
  count <- t(count)
  sub_design <- as.data.frame(phyloseq::sample_data(ps))
  dim(sub_design)
  sub_design$SampleType = as.character(sub_design$Group)
  sub_design$SampleType <- as.factor(sub_design$Group)
  # create DGE list
  d = edgeR::DGEList(counts=count, group=sub_design$SampleType)
  d$samples
  d = edgeR::calcNormFactors(d,method=method)#默认为TMM标准化

  # Building experiment matrix
  design.mat = model.matrix(~ 0 + d$samples$group)
  colnames(design.mat)=levels(sub_design$SampleType)
  d2 = edgeR::estimateGLMCommonDisp(d, design.mat)
  d2 = edgeR::estimateGLMTagwiseDisp(d2, design.mat)
  fit = edgeR::glmFit(d2, design.mat)




  #------------根据分组提取需要的差异结果#------------
  plot_list <- list()
  for (i in 1:dim(aaa)[2]) {
    # i = 1
    Desep_group = aaa[,i]
    print( Desep_group)


    # head(design)
    # 设置比较组写在前面的分组为enrich表明第一个分组含量高

    group = paste(Desep_group[1],Desep_group[2],sep = "-")
    group
    BvsA <- limma::makeContrasts(contrasts =  group,levels=design.mat)#注意是以GF1为对照做的比较
    # 组间比较,统计Fold change, Pvalue
    lrt = edgeR::glmLRT(fit,contrast=BvsA)

    # FDR检验，控制假阳性率小于5%
    de_lrt = edgeR::decideTestsDGE(lrt, adjust.method="fdr", p.value=pvalue,lfc=lfc)#lfc=0这个是默认值
    summary(de_lrt)
    # 导出计算结果
    x=lrt$table
    x$sig=de_lrt
    head(x)
    #------差异结果符合otu表格的顺序
    row.names(count)[1:6]

    x <- cbind(x, padj = p.adjust(x$PValue, method = "fdr"))
    enriched = row.names(subset(x,sig==1))
    depleted = row.names(subset(x,sig==-1))

    x$level = as.factor(ifelse(as.vector(x$sig) ==1, "enriched",ifelse(as.vector(x$sig)==-1, "depleted","nosig")))
    x = data.frame(row.names = row.names(x),logFC = x$logFC,level = x$level,p = x$PValue)
    head(x)
    # colnames(x) = paste(group,colnames(x),sep = "")



    # x = res
    # head(x)
    #------差异结果符合otu表格的顺序
    # x = data.frame(row.names = row.names(x),logFC = x$log2FoldChange,level = x$level,p = x$pvalue)
    x1 = x %>%
      dplyr::filter(level %in% c("enriched","depleted","nosig") )
    head(x1)
    x1$Genus = row.names(x1)
    # x$level = factor(x$level,levels = c("enriched","depleted","nosig"))
    if (nrow(x1)<= 1) {

    }
    x2 <- x1 %>%
      dplyr::mutate(ord = logFC^2) %>%
      dplyr::filter(level != "nosig") %>%
      dplyr::arrange(desc(ord)) %>%
      head(n = 5)
    x3 <- x1 %>%
      dplyr::mutate(ord = logFC^2) %>%
      dplyr::filter(level != "nosig") %>%
      dplyr::arrange(desc(ord))
    #
    # file = paste(path,"/",group,j,"_","Edger_Volcano_difference.csv",sep = "")
    # write.csv(x3,file,quote = F)
    head(x2)

    p <- ggplot(x1,aes(x =logFC ,y = -log2(p), colour=level)) +
      geom_point() +
      geom_hline(yintercept=-log10(0.2),
                 linetype=4,
                 color = 'black',
                 size = 0.5) +
      geom_vline(xintercept=c(-1,1),
                 linetype=3,
                 color = 'black',
                 size = 0.5) +
      ggrepel::geom_text_repel(data=x2, aes(x =logFC ,y = -log2(p), label=Genus), size=1) +
      scale_color_manual(values = c('blue2','red2', 'gray30')) +
      ggtitle(group) + theme_bw()

    p

    # file = paste(path,"/",group,j,"_","Edger_Volcano.pdf",sep = "")
    # ggsave(file,p,width = 8,height = 6)
    #
    # file = paste(path,"/",group,j,"_","Edger_Volcano.png",sep = "")
    # ggsave(file,p,width = 8,height = 6)
    #
    plot_list[[i]] <- p
    colnames(x) = paste(group,colnames(x),sep = "")



    if (i ==1) {
      table =x
    }
    if (i != 1) {
      table = cbind(table,x)
    }
  }

  x = table

  ###########添加物种丰度#----------
  # dim(count)
  # str(count)
  count = as.matrix(count)
  norm = t(t(count)/colSums(count)) #* 100 # normalization to total 100
  dim(norm)
  norm1 = norm %>%
    t() %>% as.data.frame()
  # head(norm1)
  #数据分组计算平均值
  library("tidyverse")
  head(norm1)

  iris.split <- split(norm1,as.factor(sub_design$SampleType))
  iris.apply <- lapply(iris.split,function(x)colMeans(x))
  # 组合结果
  iris.combine <- do.call(rbind,iris.apply)
  norm2= t(iris.combine)

  #head(norm)
  str(norm2)
  norm2 = as.data.frame(norm2)
  # dim(x)
  head(norm2)
  x = cbind(x,norm2)
  head(x)
  #在加入这个文件taxonomy时，去除后面两列不相干的列
  # 读取taxonomy，并添加各列名称

  if (!is.null(ps@tax_table)) {
    taxonomy = as.data.frame(ggClusterNet::vegan_tax(ps))
    head(taxonomy)
    # taxonomy <- as.data.frame(tax_table(ps1))

    #发现这个注释文件并不适用于直接作图。
    #采用excel将其分列处理，并且删去最后一列，才可以运行
    if (length(colnames(taxonomy)) == 6) {
      colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus")
    }else if (length(colnames(taxonomy)) == 7) {
      colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species")
    }else if (length(colnames(taxonomy)) == 8) {
      colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species","rep")
    }
    # colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus")

    # Taxonomy排序，并筛选OTU表中存在的
    library(dplyr)
    taxonomy$id=rownames(taxonomy)
    # head(taxonomy)
    tax = taxonomy[row.names(x),]
    x = x[rownames(tax), ] # reorder according to tax

    if (length(colnames(taxonomy)) == 7) {
      x = x[rownames(tax), ] # reorder according to tax
      x$phylum = gsub("","",tax$phylum,perl=TRUE)
      x$class = gsub("","",tax$class,perl=TRUE)
      x$order = gsub("","",tax$order,perl=TRUE)
      x$family = gsub("","",tax$family,perl=TRUE)
      x$genus = gsub("","",tax$genus,perl=TRUE)
      # x$species = gsub("","",tax$species,perl=TRUE)
    }else if (length(colnames(taxonomy)) == 8) {
      x$phylum = gsub("","",tax$phylum,perl=TRUE)
      x$class = gsub("","",tax$class,perl=TRUE)
      x$order = gsub("","",tax$order,perl=TRUE)
      x$family = gsub("","",tax$family,perl=TRUE)
      x$genus = gsub("","",tax$genus,perl=TRUE)
      x$species = gsub("","",tax$species,perl=TRUE)
    }else if (length(colnames(taxonomy)) == 9) {
      x = x[rownames(tax), ] # reorder according to tax
      x$phylum = gsub("","",tax$phylum,perl=TRUE)
      x$class = gsub("","",tax$class,perl=TRUE)
      x$order = gsub("","",tax$order,perl=TRUE)
      x$family = gsub("","",tax$family,perl=TRUE)
      x$genus = gsub("","",tax$genus,perl=TRUE)
      x$species = gsub("","",tax$species,perl=TRUE)


    }

  } else {
    x = cbind(x,tax)
  }

  return(list(plot_list,x))

}

