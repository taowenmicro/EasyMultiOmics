#' @title Perform Piecewise Structural Equation Modeling (SEM) for Omics Data
#' @description
#' This function performs Piecewise Structural Equation Modeling (SEM) on multiple omics datasets (e.g., transcriptomics, metabolomics, microbiome, metagenome).
#' The function uses Principal Component Analysis (PCA) to extract the first principal component from each omics dataset, filters low abundance taxa,
#' and then builds a piecewise SEM model to examine the relationships between different omics layers.
#'
#' The function takes multiple omics datasets (transcriptome, metabolome, microbiome, and metagenome), intersects their samples,
#' filters low abundance features, and performs PCA. Then, a piecewise SEM model is constructed to understand the relationships between
#' these datasets, such as the influence of transcriptomics on metabolomics, metabolomics on microbiomes, and the combined effects of
#' metabolomics and microbiomes on metagenomes. The function returns the R², AIC, and BIC values of the model.
#' @return
#' A list containing:
#' - `r_squared`: The R² value of the model.
#' - `aic`: The Akaike Information Criterion (AIC) value of the model.
#' - `bic`: The Bayesian Information Criterion (BIC) value of the model.
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
#' res <- piecewiseSEM.omics(ps.trans=ps.trans,ps.ms=ps.ms,ps.micro=ps.micro,
#' ps.meta=ps.kegg,filter=0.05)
#' res
#' }
#'
piecewiseSEM.omics <-function(ps.trans=ps.trans,ps.ms=ps.ms,ps.micro=ps.micro,
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
  pca_data <- scale(pca_data) %>%as.data.frame()
  ##3 #开始piecewiseSEM建模----
  mod1 <- piecewiseSEM::psem(
    lm(Metabolome ~ Transcriptome, data = pca_data),
    lm(Microbiome ~ Metabolome, data = pca_data),
    lm(Metagenome ~ Metabolome + Microbiome, data = pca_data)
  )# print(summary(sem1,standardized = T, rsq = T))
  print(plot(mod1))
  # 提取模型摘要
  model_summary <- summary(mod1)

  # 提取 R²、AIC、BIC 等统计值
  r_squared <- model_summary$R2
  aic <- model_summary$AIC
  bic <- model_summary$BIC

  # 将统计值格式化为字符串
  # fit_text <- paste(
  #   "R-squared = ", round(r_squared, 2), "\n",
  #   "AIC = ", round(aic, 2), "\n",
  #   "BIC = ", round(bic, 2)
  # )

  return(list(r_squared,aic,bic))
}
