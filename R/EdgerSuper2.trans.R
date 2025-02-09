#' @title Differential Abundance Analysis Using EdgeR for Groups
#'
#' @description
#' The `EdgerSuper2.micro` function performs differential abundance analysis between microbial groups using the `EdgeR` package.
#' It calculates log fold changes, p-values, and adjusted p-values (FDR) to identify significantly enriched or depleted taxa
#' between group comparisons. The function works on data stored in a `phyloseq` object or manually provided OTU, taxonomic, and metadata tables.
#'
#' @param otu A data frame containing OTU counts. Optional if `ps` is provided.
#' @param tax A data frame containing taxonomic annotations. Optional if `ps` is provided.
#' @param map A data frame containing sample metadata. Optional if `ps` is provided.
#' @param tree A phylogenetic tree. Optional if `ps` is provided.
#' @param ps A `phyloseq` object containing microbiome data. If provided, overrides `otu`, `tax`, `map`, and `tree`.
#' @param group A string specifying the grouping variable in the sample metadata. Default is `"Group"`.
#' @param pvalue A numeric value specifying the significance threshold for adjusted p-values (FDR). Default is `0.05`.
#' @param lfc A numeric value specifying the log2 fold change cutoff for significance. Default is `0`.
#' @param artGroup A custom matrix specifying pairwise group comparisons. Optional.
#' @param method A string specifying the normalization method for EdgeR. Default is `"TMM"`.
#' @param j A string or integer specifying the taxonomic rank to perform the analysis on.
#' Can be a numeric rank (1-7), a taxonomic name (e.g., `"Phylum"`), or `"OTU"`. Default is `2`.
#'
#' @return
#' A data frame containing the differential abundance results:
#' \describe{
#'   \item{logFC}{The log2 fold change of taxa between the compared groups.}
#'   \item{p}{The raw p-values for each taxon.}
#'   \item{padj}{The adjusted p-values (FDR).}
#'   \item{level}{The classification of taxa as `"enriched"`, `"depleted"`, or `"nosig"` based on the logFC and p-value thresholds.}
#'   \item{group}{The name of the group comparison (e.g., `"Group1-Group2"`).}
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Prepares OTU, taxonomic, and metadata tables from the `phyloseq` object or individual inputs.
#'   \item Aggregates data to the specified taxonomic rank (if applicable).
#'   \item Normalizes OTU counts using the specified EdgeR normalization method (`"TMM"` by default).
#'   \item Constructs a model matrix and fits a generalized linear model (GLM) to the data.
#'   \item Performs pairwise group comparisons using contrasts in the EdgeR framework.
#'   \item Calculates log fold changes, raw p-values, and FDR-adjusted p-values for each taxon.
#'   \item Labels taxa as `"enriched"`, `"depleted"`, or `"nosig"` based on thresholds.
#'   \item Outputs a combined data frame with differential abundance results and group comparisons.
#' }
#'
#' The results can be visualized using additional plotting tools, such as volcano plots, which highlight the significant taxa.
#'
#' @examples
#' \dontrun{
#' res =  EdgerSuper2.trans (ps = ps.trans,group  = "Group",artGroup =NULL, j = "gene")
#' head(res)
#'
#' }
#'
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#'
#' @export
EdgerSuper2.trans= function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,
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

    # file = paste(path,"/",group,j,"_","Edger_Volcano_Top5.csv",sep = "")
    # write.csv(x2,file,quote = F)
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

    # p
    #
    # file = paste(path,"/",group,j,"_","Edger_Volcano.pdf",sep = "")
    # ggsave(file,p,width = 8,height = 6)
    #
    # file = paste(path,"/",group,j,"_","Edger_Volcano.png",sep = "")
    # ggsave(file,p,width = 8,height = 6)


    colnames(x) = paste(colnames(x),sep = "")
    x$group = group


    if (i ==1) {
      table =x
    }
    if (i != 1) {
      table = rbind(table,x)
    }
  }
  x = table
  return(x)

}

