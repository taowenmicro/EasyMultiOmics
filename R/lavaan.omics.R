#' @title Multi-omics structural equation model construction--Lavaan
#' @description
#' This function is used to conduct Structural Equation Modeling (SEM) analysis by lavaan,
#' the relationship between Transcriptome,Metabolome,Microbiome and Metagenome is explored by extracting the results of the first axis of PCA analysis.
#' @param ps.trans Transcriptome data for analysis, a phyloseq format file used as an alternative for the input containing transcriptome functional composition table.
#' @param ps.ms Metabolome data for analysis, a phyloseq format file used as an alternative for the input containing metabolites composition table.
#' @param ps.micro Microbiome data for analysis, a phyloseq format file used as an alternative for the input containing microbial composition table.
#' @param ps.meta Metagenome functional data for analysis, a phyloseq format file used as an alternative for the input containing metagenome functional composition table.
#' @param filter Threshold for filtering low abundance data,default is 0.001.
#'
#' @return Returns a data frame containing model fit statistics.
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njau.edu.cn}.
#' @examples
#' \dontrun{
#' library(vegan)
#' library(phyloseq)
#' library(dplyr)
#' library(ggClusterNet)
#' library(lavaan)
#' library(semPlot)
#' library(EasyMultiOmics)
#' data(ps.trans)
#' data(ps.ms)
#' data(ps.micro)
#' data(ps.kegg)
#' res <- lavaan.omics(ps.trans=ps.trans,ps.ms=ps.ms,ps.micro=ps.micro,
#' ps.meta=ps.kegg,filter=0.05)
#' res
#' }
lavaan.omics <-function(ps.trans=ps.trans,ps.ms=ps.ms,ps.micro=ps.micro,
                        ps.meta=ps.kegg,filter=0.001){
  ##1 共同存在样本----
  data1 <- data.frame(otu_table(ps.trans))
  data2 <- data.frame(otu_table(ps.ms))
  data3 <- data.frame(otu_table(ps.micro))
  data4 <- data.frame(otu_table(ps.kegg))
  union_sample <- base::Reduce(intersect, list(colnames(data1),colnames(data2),colnames(data3),colnames(data4)))
  print(paste("union_sample_number: ",length(union_sample),sep = ""))
  ##2 去除低丰度数据,提取pca1轴数据----
  ps_list <- list(ps.trans,ps.ms,ps.micro,ps.kegg)
  pca_data <- data.frame()
  #ps <- ps.kegg
  for(ps in ps_list){
    tem = ps %>% phyloseq::transform_sample_counts(function(x) x/sum(x)) %>%
      phyloseq::filter_taxa(function(x) sum(x) > filter, TRUE) %>%
      ggClusterNet::vegan_otu()
    otu_table <- tem[union_sample,]
    otu.pca = stats::prcomp(otu_table, scale. = TRUE)
    points = data.frame(otu.pca$x[,1])
    col_name <- paste("ps_", ncol(pca_data)+ 1, sep = "")
    if (ncol(pca_data) == 0) {
      pca_data<- points
      colnames(pca_data) <- col_name
    } else {
      pca_data[[col_name]] <- points[,1]
    }
  }
  colnames(pca_data) <- c("Transcriptome","Metabolome","Microbiome","Metagenome")
  pca_data <- scale(pca_data)
  ##3 #开始lavaan建模----
  mod1 <- " Metabolome~ Transcriptome
         Microbiome~ Metabolome
         Metagenome~Metabolome+Microbiome"

  sem1 <- lavaan::sem(mod1, data = pca_data)
  # print(summary(sem1,standardized = T, rsq = T))
  model_measures <- data.frame(fitMeasures=lavaan::fitMeasures(sem1,
                                                               c("chisq","df","pvalue","gfi","cfi","rmr","srmr","rmsea")))

  semPlot::semPaths(sem1,
                    what = "std",
                    layout = "tree2", #图的格式， tree  circle  spring  tree2  circle2
                    fade = F, #褪色，按照相关度褪色
                    residuals = F ,#残差/方差要不要体现在图中，可T和F
                    nCharNodes = 0,
                    sizeMan= 12,
                    edge.label.cex=0.7)
  return(model_measures)
}
