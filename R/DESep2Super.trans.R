#' @title Differential Abundance Analysis Using DESeq2
#'
#' @description
#' The `DESep2Super.metf` function performs differential abundance analysis between transcriptome functional groups using the `DESeq2` package.
#' It identifies significantly enriched or depleted taxa across multiple group comparisons by calculating log2 fold changes,
#' raw p-values, and adjusted p-values (FDR). The input data can be provided as a `phyloseq` object or as separate transcriptome functional composition, taxonomic, and metadata tables.
#'
#' @param otu A data frame containing gene counts. Optional if `ps` is provided.
#' @param map A data frame containing sample metadata. Optional if `ps` is provided.
#' @param tree A phylogenetic tree. Optional if `ps` is provided.
#' @param tax  transcriptome functional classification table.Optional if `ps` is provided.
#' @param ps A phyloseq format file used as an alternative for the input containing transcriptome functional composition table, tax, and sample metadata.
#' @param j A string or integer specifying the taxonomic rank to perform the analysis on.
#' Can be a numeric rank (1-7), a taxonomic name (e.g., `"Phylum"`), or `"OTU"`. Default is `"meta"`.
#' @param group A string specifying the grouping variable in the sample metadata. Default is `"Group"`.
#' @param pvalue A numeric value specifying the significance threshold for adjusted p-values (FDR). Default is `0.05`.
#' @param artGroup A custom matrix specifying pairwise group comparisons. Optional.
#' @return
#' A list containing:
#' \describe{
#'   \item{A list of volcano plots for each pairwise comparison.}
#'   \item{A data frame containing differential abundance results, including log2 fold changes, p-values, adjusted p-values,
#'   and group-specific normalized abundances.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Prepares OTU, taxonomic, and metadata tables from the `phyloseq` object or individual inputs.
#'   \item Aggregates data to the specified taxonomic rank (if applicable).
#'   \item Constructs a `DESeq2` dataset object and normalizes OTU counts.
#'   \item Performs pairwise group comparisons using contrasts in the DESeq2 framework.
#'   \item Calculates log2 fold changes, raw p-values, and FDR-adjusted p-values for each taxon.
#'   \item Labels taxa as `"enriched"`, `"depleted"`, or `"nosig"` based on thresholds.
#'   \item Outputs a combined data frame with differential abundance results, normalized abundances, and taxonomic annotations.
#' }
#'
#' The results include volcano plots that visualize the differential abundance results for each group comparison.
#' Significant taxa are highlighted, and their log2 fold changes and p-values are displayed.
#'
#' @examples
#' \dontrun{
#' data(ps.trans)
#' ps =ps.trans %>% filter_OTU_ps(Top = 1000)
#' res = DESep2Super.trans(ps = ps,group  = "Group",artGroup = NULL)
#' p15.1 = res[[1]][1]
#' p15.1
#' p15.2 = res[[1]][2]
#' p15.2
#' p15.3 = res[[1]][3]
#' p15.3
#' dat = res[[2]]
#' head(dat)
#' }
#'
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#'
#' @export
DESep2Super.trans= function(otu = NULL,
                            tax = NULL,
                            map = NULL,
                            tree = NULL ,
                            ps = NULL,
                            j = "meta",
                            group  = "Group",
                            pvalue = 0.05,
                            artGroup = NULL
){


  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)
  # ps = ps %>%
  #   ggClusterNet::tax_glom_wt(ranks = j)
  if (j %in% c("OTU","gene","meta","trans")) {
    ps = ps
  } else if (j %in% c(1:7)) {
    ps = ps %>%
      ggClusterNet::tax_glom_wt(ranks = j)
  } else if (j %in% c("Kingdom","Phylum","Class","Order","Family","Genus","Species")){

  } else {
    ps = ps
    print("unknown j, checked please")
  }


  Desep_group <- ps %>%
    phyloseq::sample_data() %>%
    .$Group %>%
    as.factor() %>%
    levels() %>%
    as.character()

  if ( is.null(artGroup)) {
    #--构造两两组合#-----
    aaa = combn(Desep_group,2)
    # sub_design <- as.data.frame(sample_data(ps))
  } else if (!is.null(artGroup)) {
    aaa  = as.matrix(b )
  }

  count <-  ps %>%
    ggClusterNet::vegan_otu() %>% round(0) %>%
    t()


  map = ps %>%
    phyloseq::sample_data() %>%
    as.tibble() %>%
    as.data.frame()


  dds <- DESeq2::DESeqDataSetFromMatrix(countData = count,
                                        colData = map,
                                        design = ~ Group)

  dds2 <- DESeq2::DESeq(dds)  ##第二步,标准化
  # resultsNames(dds2)


  #------------根据分组提取需要的差异结果#------------
  plot_list <- list()
  for (i in 1:dim(aaa)[2]) {
    # i = 1
    Desep_group = aaa[,i]
    # print( Desep_group)


    # head(design)
    # 设置比较组写在前面的分组为enrich表明第一个分组含量高

    group = paste(Desep_group[1],Desep_group[2],sep = "-")
    group

    # 将结果用results()函数来获取，赋值给res变量
    res <-  DESeq2::results(dds2, contrast=c("Group",Desep_group ),alpha=0.05)
    # 导出计算结果

    x = res
    head(x)
    #------差异结果符合otu表格的顺序

    x$level = as.factor(ifelse(as.vector(x$padj) < 0.05 & x$log2FoldChange > 0, "enriched",
                               ifelse(as.vector(x$padj) < 0.05 &x$log2FoldChange < 0, "depleted","nosig")))


    x = data.frame(row.names = row.names(x),logFC = x$log2FoldChange,level = x$level,p = x$pvalue)
    x1 = x %>%
      filter(level %in% c("enriched","depleted","nosig") )
    head(x1)
    x1$Genus = row.names(x1)
    # x$level = factor(x$level,levels = c("enriched","depleted","nosig"))

    x2 <- x1 %>%
      dplyr::mutate(ord = logFC^2) %>%
      dplyr::filter(level != "nosig") %>%
      dplyr::arrange(desc(ord)) %>%
      head(n = 5)

    x3 <- x1 %>%
      dplyr::mutate(ord = logFC^2) %>%
      dplyr::filter(level != "nosig") %>%
      dplyr::arrange(desc(ord))

    # file = paste(path,"/",group,j,"_","DESep2_Volcano_diffference.csv",sep = "")
    # write.csv(x3,file,quote = F)
    head(x2)
    head(x1)
    x1$p[is.na(x1$p)] = 1
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
    plot_list[[i]] <- p
    # file = paste(path,"/",group,j,"_","DESep2_Volcano.pdf",sep = "")
    # ggsave(file,p,width = 8,height = 6)
    #
    # file = paste(path,"/",group,j,"_","DESep2_Volcano.png",sep = "")
    # ggsave(file,p,width = 8,height = 6)

    colnames(x) = paste(group,colnames(x),sep = "")





    if (i ==1) {
      table =x
    }
    if (i != 1) {
      table = cbind(table,x)
    }
  }

  # x = table

  # head(table)
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
  # library("tidyverse")
  head(norm1)

  iris.split <- split(norm1,as.factor(map$Group))
  iris.apply <- lapply(iris.split,function(x)colMeans(x))
  # 组合结果
  iris.combine <- do.call(rbind,iris.apply)
  norm2= t(iris.combine)

  #head(norm)
  str(norm2)
  norm2 = as.data.frame(norm2)

  head(norm2)
  x = cbind(table,norm2)
  head(x)
  #在加入这个文件taxonomy时，去除后面两列不相干的列
  # 读取taxonomy，并添加各列名称
  #
  # if (!is.null(ps@tax_table)) {
  #   taxonomy = as.data.frame(ggClusterNet::vegan_tax(ps))
  #   head(taxonomy)
  #   # taxonomy <- as.data.frame(tax_table(ps1))
  #
  #   #发现这个注释文件并不适用于直接作图。
  #   #采用excel将其分列处理，并且删去最后一列，才可以运行
  #   if (length(colnames(taxonomy)) == 6) {
  #     colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus")
  #   }else if (length(colnames(taxonomy)) == 7) {
  #     colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species")
  #   }else if (length(colnames(taxonomy)) == 8) {
  #     colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species","rep")
  #   }
  #   # colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus")
  #
  #   # Taxonomy排序，并筛选OTU表中存在的
  #   library(dplyr)
  #   taxonomy$id=rownames(taxonomy)
  #   # head(taxonomy)
  #   tax = taxonomy[row.names(x),]
  #   x = x[rownames(tax), ] # reorder according to tax
  #
  #   if (length(colnames(taxonomy)) == 7) {
  #     x = x[rownames(tax), ] # reorder according to tax
  #     x$phylum = gsub("","",tax$phylum,perl=TRUE)
  #     x$class = gsub("","",tax$class,perl=TRUE)
  #     x$order = gsub("","",tax$order,perl=TRUE)
  #     x$family = gsub("","",tax$family,perl=TRUE)
  #     x$genus = gsub("","",tax$genus,perl=TRUE)
  #     # x$species = gsub("","",tax$species,perl=TRUE)
  #   }else if (length(colnames(taxonomy)) == 8) {
  #     x$phylum = gsub("","",tax$phylum,perl=TRUE)
  #     x$class = gsub("","",tax$class,perl=TRUE)
  #     x$order = gsub("","",tax$order,perl=TRUE)
  #     x$family = gsub("","",tax$family,perl=TRUE)
  #     x$genus = gsub("","",tax$genus,perl=TRUE)
  #     x$species = gsub("","",tax$species,perl=TRUE)
  #   }else if (length(colnames(taxonomy)) == 9) {
  #     x = x[rownames(tax), ] # reorder according to tax
  #     x$phylum = gsub("","",tax$phylum,perl=TRUE)
  #     x$class = gsub("","",tax$class,perl=TRUE)
  #     x$order = gsub("","",tax$order,perl=TRUE)
  #     x$family = gsub("","",tax$family,perl=TRUE)
  #     x$genus = gsub("","",tax$genus,perl=TRUE)
  #     x$species = gsub("","",tax$species,perl=TRUE)
  #
  #
  #   }
  #
  # } else {
  #   x = x
  # }
  #
  return(list(plot_list,x))

}



