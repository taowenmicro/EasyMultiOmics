#
# diffpath = paste(otupath,"/diff_Manhattan/",sep = "")
# dir.create(diffpath)
#
# edge_Manhattan (
#   ps = ps,
#   pvalue = 0.05,
#   lfc = 0,
#   diffpath = diffpath
# )


#--曼哈顿图
# library("ggplot2")
# library("vegan")
# library("RColorBrewer")
# library(edgeR)

edge_Manhattan.micro <- function(
  ps = ps,
  pvalue = 0.05,
  lfc = 0
){

  count = ps %>%
    ggClusterNet::vegan_otu() %>%
    t()
  # create DGE list

  d = edgeR::DGEList(counts=count, group= phyloseq::sample_data(ps)$Group)
  d = edgeR::calcNormFactors(d)
  design.mat = model.matrix(~ 0 + d$samples$group)
  colnames(design.mat)=levels(as.factor(phyloseq::sample_data(ps)$Group))
  d2 = edgeR::estimateGLMCommonDisp(d, design.mat)
  d2 = edgeR::estimateGLMTagwiseDisp(d2, design.mat)
  fit = edgeR::glmFit(d2, design.mat)

  Desep_group <- as.character(levels(as.factor(phyloseq::sample_data(ps)$Group)))
  Desep_group

  aaa = combn(Desep_group,2)
  plot_list <- list()
  for (i in 1:dim(aaa)[2]) {
    # i = 1
    Desep_group = aaa[,i]
    print( Desep_group)
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


    x$sig=de_lrt
    head(x)
    #------差异结果符合otu表格的顺序
    row.names(count)[1:6]

    x <- cbind(x, padj = p.adjust(x$PValue, method = "fdr"))
    enriched = row.names(subset(x,sig==1))
    depleted = row.names(subset(x,sig==-1))
    x$level = as.factor(ifelse(as.vector(x$sig) ==1, "enriched",ifelse(as.vector(x$sig)==-1, "depleted","nosig")))
    x$otu = rownames(x)
    x$neglogp = -log(x$PValue)

    tax = ps %>%
      ggClusterNet::vegan_tax() %>%
      as.data.frame()
    head(tax)
    x = cbind(x,tax)
    head(x)

    top_phylum=c("Bacteroidetes","Firmicutes","Planctomycetes","Proteobacteria","Verrucomicrobia")
    x[!(x$Phylum %in% top_phylum),]$Phylum = "Low Abundance" # no level can get value
    x$otu = factor(x$otu, levels=x$otu)   # set x order
    x$level = factor(x$level, levels=c("enriched","depleted","nosig"))
    levels(x$Phylum)=c(top_phylum,"Low Abundance")
    x = x %>% arrange(Phylum)
    head(x)
    x$otu = factor(x$otu, levels=x$otu)




    # if (tem != 0) {
    #   tem = x[x$neglogp>15,] %>% nrow()
    # }

    FDR = min(x$neglogp[x$level=="depleted"])
    p = ggplot(x, aes(x=otu, y=neglogp, color=Phylum, size=logCPM, shape=level)) +
      geom_point(alpha=.7) +
      geom_hline(yintercept=FDR, linetype=2, color="lightgrey") +
      scale_shape_manual(values=c(17, 25, 20))+
      scale_size(breaks=c(5, 10, 15)) +
      labs(x="OTU", y="-loge(P)") +
      theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),legend.position="top") +
      scale_color_manual(values = c(RColorBrewer::brewer.pal(9,"Set1")))+
      ggtitle(group)
    p

    # filename = paste(diffpath,"/",paste(Desep_group[1],Desep_group[2],sep = "_"),"Manhattan_plot.pdf",sep = "")
    # ggsave(filename,p,width = 16,height = 6)
    # filename = paste(diffpath,"/",paste(Desep_group[1],Desep_group[2],sep = "_"),"Manhattan_plot.png",sep = "")
    # ggsave(filename,p,width = 16,height = 6,dpi = 72)
    plot_list[[i]] <- p
  }
  return(plot_list)
}











