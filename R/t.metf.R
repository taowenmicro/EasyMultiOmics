#' @title Differential abundance analysis using t test for metagenome functional data
#' @description
#' This function performs differential abundance analysis for metagenome functional data using the `t.test` function.
#' It identifies differentially abundant taxa betwee two groups in a `phyloseq` object or selected comparison groups,
#' visualizes the results with volcano plots, and outputs realtive abundance data with test results.
#' @param ps A phyloseq format file used as an alternative for the input containing metagenome
#' functional composition table, tax, and sample metadata.
#' @param artGroup A matrix specifying custom group comparisons. Optional,
#' the standard format is two rows,Each column represents the two groups selected for comparison,
#' with the first row defaulting to the reference control group.
#' @param group A string specifying the grouping variable in the sample metadata. Default is `"Group"`.
#' @param pvalue A numeric value specifying the significance threshold for p value. Default is `0.05`.
#' @param FCvalue A numeric value specifying the significance threshold for fold change value. Default is `2`.
#' @param padjust A logical value that is used to determine whether to perform p-value correction. Default is `TRUE`.
#' @param padjust_method Select the method for p-value correction("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#' "fdr"),when padjust is `TRUE`. Default is `BH`.
#' @param log2FC_threshold A numeric value specifying the plot threshold for
#' log2fold_change_value. Default is `10`.
#' @param ncol Number of columns for arranging plots.Default is `3`.
#' @param nrow Number of rows for arranging plots.Default is `1`.
#' @param outpath Path to save the analysis results and plots.
#' @author
#' Tao Wen \email{taowen@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njau.edu.cn}
#' @return A list containing the following components:
#' \item{Volcano_all_plot}{A combined volcano plots of all comparison groups.}
#' \item{result_all}{Data frame containing all comparison groups analysis result.}
#' \item{meta}{Data frame containing the union and intersection of differentially expressed genes across all pairwise comparison groups.}
#' \item{plot_list}{A list of individual volcano plots.}
#' @export
#' @examples
#' \dontrun{
#' ps=ps.kegg %>% filter_OTU_ps(Top = 1000)
#' res = t_metf(ps = ps,group  = "Group",artGroup =NULL,group  = "Group",
#'              pvalue=0.05,FCvalue=2,padjust=TRUE,padjust_method="BH",
#'              log2FC_threshold=10,ncol=3,nrow=1,outpath=NULL)
#' p17=res[[1]]
#' p17
#' data=res[[2]]
#' head(data)
#' inter_union=res[[3]]
#' head(inter_union)
#' p17.1 = res [[4]][1]
#' p17.1
#' p17.2 = res [[4]][2]
#' p17.2
#' p17.3 = res [[4]][3]
#' p17.3
#' }

t_metf = function(
    ps = ps,
    artGroup = NULL,
    group  = "Group",
    pvalue=0.05,
    FCvalue=2,
    padjust=TRUE,
    padjust_method="BH",
    log2FC_threshold=10,
    ncol=3,
    nrow=1,
    outpath=NULL
){
#Compare group settings
sub_design <- as.data.frame(phyloseq::sample_data(ps))
Desep_group <- as.character(levels(as.factor(sub_design[[group]])))
if ( is.null(artGroup)) {
  aaa = combn(Desep_group,2)
} else if (!is.null(artGroup)) {
  aaa  = as.matrix(artGroup)
}
#Gene abundance calculations
count = t(ggClusterNet::vegan_otu(ps))%>%as.matrix()
norm = t(t(count)/colSums(count))#*100 # normalization to total 100
norm1 = norm %>%t() %>% as.data.frame()
##Group gene abundance calculations
iris.split <- split(norm1,as.factor(sub_design[[group]]))
iris.apply <- lapply(iris.split,function(x) colMeans(x))
iris.combine <- do.call(rbind,iris.apply)
norm2= as.data.frame(t(iris.combine))
norm2$ID=rownames(norm2)
plot_list <- list()
for (i in 1:dim(aaa)[2]) {
  Desep_group = aaa[,i]
  group1 = paste(Desep_group[2],Desep_group[1],sep = "__")
  ##5.3构建对应ps文件，转换为相对丰度----
  ps11 <-ps %>%subset_samples.wt(group,c(Desep_group[2] ,Desep_group[1]))%>%
    phyloseq::transform_sample_counts(function(x) x / sum(x)*10000 )
  ##5.4构建计算数据
  data <- t(data.frame(otu_table(ps11)))
  tax1 <- norm2[colnames(data),]
  sample_group <- sample_data(ps11)
  sample_group <- sample_group[rownames(data),]
  data <- cbind(sample_group[,group],data)
  ##5.5计算差异倍数----
  log2FC <-log2((tax1[,Desep_group[2]]+1/mean(colSums(count))) /(tax1[,Desep_group[1]]+1/mean(colSums(count))))
  #差异检验
  result=data.frame(ID=colnames(data)[2:ncol(data)],
                    p=sapply(data[,2:ncol(data)],function(x){t.test(x~data$Group)[["p.value"]]}),
                    logFC=log2FC)
  #padjust
  result <- cbind(result,padj=p.adjust(result$p, method =padjust_method))
  if(padjust==TRUE){
  result$level = as.factor(ifelse(as.vector(result$logFC)>log2(FCvalue)&as.vector(result$padj)<pvalue, "enriched",
                          ifelse(as.vector(result$logFC)<(-log2(FCvalue))&as.vector(result$padj)<pvalue, "depleted",
                                 "nosig")))}else{
  result$level= as.factor(ifelse(as.vector(result$logFC)>log2(FCvalue)&as.vector(result$p)<pvalue, "enriched",
                              ifelse(as.vector(result$logFC)<(-log2(FCvalue))&as.vector(result$p)<pvalue, "depleted",
                                     "nosig")))}

  result1 = result %>%
    dplyr::filter(level%in% c("enriched","depleted","nosig") )
  head(result1)
  result1$Genus = row.names(result1)
  result1_output <- result1[order(result1$level),]
  result1_output <- cbind(tax1[rownames(result1_output),c(Desep_group[2],Desep_group[1])],result1_output)

  result2 <- result1 %>%
    dplyr::mutate(ord = logFC^2) %>%
    dplyr::filter(level!= "nosig") %>%
    dplyr::arrange(desc(ord)) %>%
    head(n = 5)
  result3 <- result1 %>%
    dplyr::mutate(ord = logFC^2) %>%
    dplyr::filter(level!= "nosig") %>%
    dplyr::arrange(desc(ord))
  aaaaaaa <- rownames(result3)
  #可视化
  ld <- paste("depleted",nrow(result1[result1$level=="depleted",]),sep = " ")
  lr <- paste("enriched",nrow(result1[result1$level=="enriched",]),sep = " ")
  ln <- paste("nosig",nrow(result1[result1$level=="nosig",]),sep = " ")
  #filter
  result1$logFC <-  ifelse(abs(result1$logFC)>log2FC_threshold&result1$logFC<0,-log2FC_threshold,
                           ifelse(abs(result1$logFC)>log2FC_threshold&result1$logFC>0,log2FC_threshold,result1$logFC))
  result2$logFC <- result1[rownames(result2),"logFC"]
  if(padjust==TRUE){
    p0 <-ggplot(result1,aes(x =logFC ,y = -log10(padj), colour=level))+
      geom_point(alpha=0.92,size=1.1)+
    geom_point(data=result2,aes(x =logFC ,y = -log10(padj), fill=level),shape=21,size=2,colour="black")+
      ggrepel::geom_text_repel(data=result2, aes(x =logFC ,y = -log10(padj), label=Genus),
                               show.legend=FALSE,force=5,size=3)
  }else{p0 <-ggplot(result1,aes(x =logFC ,y = -log10(p), colour=level))+
    geom_point(alpha=0.92,size=1.1) +
    geom_point(data=result2,aes(x =logFC ,y = -log10(p), fill=level),shape=21,size=2,colour="black")+
    ggrepel::geom_text_repel(data=result2, aes(x =logFC ,y = -log10(p), label=Genus),
                             show.legend=FALSE,force=5,size=3)}
  p <- p0+ggtitle(group1) +
    scale_fill_manual(values = c("enriched"="#ffad73","depleted"="#26b3ff","nosig"="grey"))+
    geom_hline(yintercept=-log10(0.05),linetype="dashed",color = 'black',size = 0.35) +
    geom_vline(xintercept=c(-log2(FCvalue),log2(FCvalue)),linetype="dashed",color = 'black',size = 0.35) +
    scale_color_manual(values=c("enriched"="#ffad73","depleted"="#26b3ff","nosig"="grey"),
                       labels=c("enriched"=lr,"depleted"=ld,"nosig"=ln)) +
    theme(panel.border=element_rect(colour= "black",fill=NA,size=0.75),
          panel.grid.minor=element_blank(),
          panel.background=element_blank(),
          plot.background=element_blank(),
          axis.title=element_text(face="bold",color="black",size=11),
          plot.title =element_text(face="bold",color="black",size=12,hjust = 0.5),
          axis.text=element_text(color="black",size=10),
          legend.background=element_blank(),
          legend.title=element_blank(),
          legend.text=element_text(face="bold",color="black",size=9.5),
           legend.spacing.x=unit(0.1,"cm"),
           legend.key.spacing.x = unit(0.1,"cm"),
          legend.key = element_blank(),
          legend.key.size =unit(0.5,"cm"),
        legend.position=c(0.99,0.99),
        legend.justification = c(0.95,0.84))+
    guides(fill="none",color = guide_legend(override.aes = list(size=2.2)))+
    labs(x="log2(FC)")#去除填充图例
  plot_list[[i]] <- p

  if (!is.null(outpath)) {
    t_testpath <- paste(outpath,group1,"/",sep="")
    fs::dir_create(t_testpath,recurse =  TRUE)
    file = paste(t_testpath,group1,"_all_genes_analysisresult.csv",sep = "")
    write.csv(result1_output,file,quote = TRUE)
    file = paste(t_testpath,group1,"_difference_genes_analysisresult.csv",sep = "")
    write.csv(result3,file,quote = TRUE)
    file = paste(t_testpath,group1,"_plotdata.csv",sep = "")
    write.csv(result1,file,quote = T)
    file = paste(t_testpath,group1,"_","Edger_Volcano.pdf",sep = "")
  ggsave(file,p,width = 11.8/2.54,height = 11.2/2.54)
  file = paste(t_testpath,group1,"_","Edger_Volcano.png",sep = "")
  ggsave(file,p,width = 11.8/2.54,height = 11.2/2.54)}

if (padjust==TRUE) {
  result_merge <- data.frame(row.names = rownames(result1),ID=result1$ID,
                             logFC=result1$logFC,p=result1$padj,level=result1$level)
}else{
  result_merge <- data.frame(row.names = rownames(result1),ID=result1$ID,
                             logFC=result1$logFC,p=result1$p,level=result1$level)
}
  colnames(result_merge)[2:4] <- paste(group1, colnames(result_merge)[2:4])
  if (i==1) {
    x <- result_merge
    bbb <- aaaaaaa
    ccc <- aaaaaaa
  }else{
    x <- merge(x,result_merge,by="ID")
    bbb <- intersect(bbb,aaaaaaa)
    ccc <- union(ccc,aaaaaaa)
  }
}
result_all <- merge(norm2,x,by="ID")
meta <- data.frame(intersect_gene=c(bbb,rep("",length(ccc)-length(bbb))),
                   union_gene=ccc)
Volcano_all_plot = ggpubr::ggarrange(plotlist = plot_list,common.legend = F,ncol = ncol,nrow = nrow)
if (!is.null(outpath)) {
  file = paste(outpath,"ttest_allgroup.csv",sep = "")
  write.csv(result_all,file,quote=T)
  file = paste(outpath,"ttest_allgroup_differentID.csv",sep = "")
  write.csv(meta,file,quote=T)
  file = paste(outpath,"ttest_allgroup_Volcano.pdf",sep = "")
  ggsave(file,Volcano_all_plot,width = 11.8*ncol/2.54,height = 11.2*nrow/2.54)
  file = paste(outpath,"ttest_allgroup_Volcano.png",sep = "")
  ggsave(file,Volcano_all_plot,width = 11.8*ncol/2.54,height = 11.2*nrow/2.54)
  }
return(list(Volcano_all_plot,result_all,meta,plot_list))
}

