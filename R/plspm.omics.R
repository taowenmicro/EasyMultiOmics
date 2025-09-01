#' @title Partial Least Squares Path Modeling (PLS-PM) for Omics Data
#' @description
#' This function performs Partial Least Squares Path Modeling (PLS-PM) on multiple omics datasets (e.g., transcriptomics, metabolomics, microbiome, and metagenome).
#' It first extracts the first principal component (PCA1) from each omics dataset, filters low-abundance taxa,
#' and then constructs a PLS-PM model to examine the relationships between different omics layers.
#' This function integrates four different omics datasets and analyzes their interdependencies using PLS-PM.
#' It follows these key steps:
#' 1. Identifies common samples across all datasets.
#' 2. Filters low-abundance taxa and extracts the first principal component (PCA1) for each dataset.
#' 3. Constructs a PLS-PM model with a predefined inner and outer model.
#' 4. Evaluates the model and returns key results, including correlation and effect size estimates.
#' @return
#' A list containing:
#' - `dat`: Summary results from the PLS-PM model.
#' - `cor`: Correlation matrix of the PLS-PM latent variables.
#' - `effect`: Estimated effects from the PLS-PM model.
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
#' res <- plspm.omics(ps.trans=ps.trans,ps.ms=ps.ms,ps.micro=ps.micro,
#' ps.meta=ps.kegg,filter=0.05)
#' res
#' }
plspm.omics <-function(ps.trans=ps.trans,ps.ms=ps.ms,ps.micro=ps.micro,
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
  ##3 #开始plspm建模----
  # 定义内部模型矩阵
  inner_model <- matrix(
    c(0, 0, 0, 0,  # Transcriptome
      1, 0, 0, 0,  # Metabolome
      0, 1, 0, 0,  # Microbiome
      0, 1, 1, 0), # Metagenome
    nrow = 4, byrow = TRUE
  )

  # 设置行名和列名
  rownames(inner_model) <- colnames(inner_model) <- c("Transcriptome",
                                                      "Metabolome", "Microbiome",
                                                      "Metagenome")

  # 定义外部模型列表
  outer_model <- list(
    Transcriptome = "Transcriptome",
    Metabolome = "Metabolome",
    Microbiome = "Microbiome",
    Metagenome = "Metagenome"
  )

  modes <- c("A", "A", "A", "A")  # 所有变量为反映型

  # 运行 PLS-PM 模型
  pls_result <-plspm:: plspm(pca_data, inner_model, outer_model, modes)
  plot(pls_result)
  # 查看模型结果
  dat =summary(pls_result)
  cor =  dat$correlations
  effect=  dat$effects

  return(list(dat,cor,effect))
}
