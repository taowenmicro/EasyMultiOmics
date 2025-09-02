#' @title Perform differential statistical analysis on metabolites data
#' @description
#' This function conducts statistical analysis on metabolites data, including
#' differential abundance analysis and relative abundance calculations.
#' @param otu A data frame containing metabolite data. If `NULL`, the `ps` object is used.
#' @param tax A data frame containing taxonomy data. If `NULL`, the `ps` object is used.
#' @param map A data frame containing sample metadata. If `NULL`, the `ps` object is used.
#' @param tree A phylogenetic tree object. If `NULL`, the `ps` object is used.
#' @param ps A phyloseq format file used as an alternative for the input containing metabolite composition table,
#' metabolite classification table, and sample metadata.
#' @param group Column name in the mapping data specifying group information. Default is "Group".
#' @param pvalue Threshold p-value for differential analysis. Default is 0.05.
#' @param artGroup User-defined group combinations.
#' @param method Statistical test method for differential analysis. Default is "wilcox" (Wilcoxon test).
#' @return Data frame containing statistical results, differential abundance, and relative abundance information.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' result1 = statSuper(ps = ps.ms,group  = "Group",artGroup = NULL,method = "wilcox")
#' head(result1)
#' result2 = statSuper(ps = ps.ms,group  = "Group",artGroup = NULL,method = "ttext")
#' head(result2)
#' @export
statSuper = function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = NULL,
         group  = "Group",pvalue = 0.05,artGroup = NULL,
         method = "wilcox" ){
  #--功能函数
  ps = inputMicro(otu,tax,map,tree,ps,group  = group)
  ps
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
  aaa
  #-将NA填充为0#
  re = phyloseq::otu_table(ps)
  #--将空缺值填充为0
  re[is.na(re)] <- 0
  otu_table(ps) <- re

  #相对丰度转化
  ps_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) );ps_rela

  #-----------------差异分析过程#-----------------------
  # otu_table = as.data.frame(vegan_otu(ps))
  count = vegan_otu(ps)
  count <- t(count)
  #数据整理形式一######
  sub_design <- as.data.frame(phyloseq::sample_data(ps))
  sub_design$ID = row.names(sub_design)

  # 转换原始数据为百分比，
  # head(norm)
  norm=as.data.frame(t(t(count)/colSums(count,na=TRUE)) * 100) # normalization to total 100
  AA=norm

  #------------根据分组提取需要的差异结果#------------
  table = NULL
  for (ii in 1:dim(aaa)[2]) {
    # ii = 1
    Desep_group = aaa[,ii]
    print( Desep_group)

    # head(design)
    #预生成2个长度与输入文件行数相同的全为0的向量，将用于存储p value和差异倍数（log2FC）
    Pvalue<-c(rep(0,nrow(AA)))
    fdr<-c(rep(0,nrow(AA)))
    log2_FC<-c(rep(0,nrow(AA)))


    # df_filter<- dplyr::filter(sub_design,SampleType %in% Desep_group)

    df_filter<- sub_design[sub_design$Group%in%Desep_group,]


    head(df_filter)

    AA = as.data.frame(AA)
    AAA = AA[as.character(df_filter$ID)]
    head(AAA)

    # 设置重复数
    rep = length(as.character(df_filter$ID))/2

    a = as.matrix(AAA)

    id1 = sub_design %>%as.tibble()  %>%dplyr::filter(Group %in% Desep_group[1]) %>% .$ID
    id2 = sub_design %>%as.tibble()  %>%dplyr::filter(Group %in% Desep_group[2]) %>% .$ID
    head(a)

    # 确保 id1 和 id2 对应的列存在于矩阵 colnames(mat)
    cols_id1 <- colnames(a)[colnames(a) %in% id1]
    cols_id2 <- colnames(a)[colnames(a) %in% id2]

    ###########开始运行脚本
    if (method == "ttext") {
      for(i in 1:nrow(a)){
        #i=1
        if(sd(a[i, cols_id1, drop = FALSE])==0&&sd(a[i, cols_id2, drop = FALSE])==0){
          Pvalue[i] <-"NA"
          log2_FC[i]<-"NA"
        }else{

        y= t.test( as.numeric(a[i, cols_id1, drop = FALSE]),as.numeric(a[i, cols_id2, drop = FALSE]))

          # y=t.test(as.numeric(a[i,(1:rep)]),as.numeric(a[i,(rep+1):(rep*2)]))

          Pvalue[i]<-y$p.value
          log2_FC[i]<-log2((mean(as.numeric(a[i, cols_id1, drop = FALSE]))+0.001)/(mean(as.numeric(a[i, cols_id2, drop = FALSE]))+0.001))
          fdr[i]=p.adjust(Pvalue[i], "BH")
        }
      }
    }


    if (method == "wilcox") {
      for(i in 1:nrow(a)){
        if(sd(a[i, cols_id1, drop = FALSE])==0&&sd(a[i, cols_id2, drop = FALSE])==0){
          Pvalue[i] <-"NA"
          log2_FC[i]<-"NA"
        }else{
          y=wilcox.test(as.numeric(a[i, cols_id1, drop = FALSE]),as.numeric(a[i, cols_id2, drop = FALSE]),exact=FALSE)
          Pvalue[i]<-y$p.value
          log2_FC[i]<-log2((mean(as.numeric(a[i, cols_id1, drop = FALSE]))+0.001)/(mean(as.numeric(a[i, cols_id2, drop = FALSE]))+0.001))
          fdr[i]=p.adjust(Pvalue[i], "BH")
        }
      }
    }


    # 在原文件后面加入log2FC，p value和FDR,共3列；
    out<-cbind(log2_FC,Pvalue,fdr)

    colnames(out) = paste(paste(Desep_group[1],Desep_group[2],sep = "_"),colnames(out),sep = "_")
    # out$tax=otu_table$compound
    head(out)
    x = as.data.frame(out)
    row.names(x) = row.names(AAA)
    head(x)
    # WT <-subset(out,fdr < 0.05 & log2_FC != "NA")
    # dim(WT)
    # head(WT)

    if (ii ==1) {
      table =x
    }
    if (ii != 1) {

      table = cbind(table,x)
    }


  }

  head(table)
  norm1 = norm %>%
    t() %>% as.data.frame()

  iris.split <- split(norm1,as.factor(sub_design$Group))
  iris.apply <- lapply(iris.split,function(x)colMeans(x))
  # 组合结果
  iris.combine <- do.call(rbind,iris.apply)
  norm2= t(iris.combine)
  str(norm2)
  norm2 = as.data.frame(norm2)

  x = cbind(AA,table)
  x = cbind(x,norm2)
  if (!is.null(ps@tax_table)) {
    taxonomy = as.data.frame(vegan_tax(ps))
    taxonomy$id= rownames(taxonomy)
    x1 = merge(x,taxonomy,by = "row.names",all.x = TRUE)
    head(x1)


  } else {
    x1 = x
  }

  return(x1)

}
