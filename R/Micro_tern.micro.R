#' @title Generate Ternary Plots for Microbial Data
#'
#' @description
#' This function generates ternary plots for visualizing microbial community data in three groups. It uses relative abundance data from a `phyloseq` object and displays OTU distribution across specified groupings.
#'
#' @param ps A `phyloseq` object containing microbiome data.
#'
#' @return
#' A list containing the following elements:
#' \describe{
#'   \item{plot}{A list of ternary plots (`ggtern` objects) for different group combinations.}
#'   \item{dataplot}{A data frame with combined OTU abundance and taxonomic information.}
#'   \item{groups}{A matrix showing all possible combinations of three groups used for ternary plots.}
#' }
#'
#' @details
#' This function performs the following steps:
#' \itemize{
#'   \item Standardizes the OTU table to relative abundances using `phyloseq::transform_sample_counts`.
#'   \item Splits the OTU data into groups based on the metadata column `Group`.
#'   \item Computes the mean abundance for each OTU within each group.
#'   \item Generates ternary plots for all possible combinations of three groups using the `ggtern` package.
#' }
#'
#' The plots visualize the relative abundances of OTUs in three groups, with points sized by their mean abundance and colored by their taxonomic classification (Phylum).
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' \dontrun{
#' ps1 = ps.16s %>% filter_OTU_ps(500)
#' res = Micro_tern.micro(ps1)
#' p15 = res[[1]]
#' dat =  res[[2]]
#' }
#'
#' @export

Micro_tern.micro = function(ps = ps
){

  ps_rela = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) );ps_rela

  otu = ggClusterNet::vegan_otu(ps_rela) %>% as.data.frame()

  #数据分组
  iris.split <- split(otu,as.factor(as.factor(phyloseq::sample_data(ps)$Group)))
  #数据分组计算平均值
  iris.apply <- lapply(iris.split,function(x)colSums(x[]))
  # 组合结果
  iris.combine <- do.call(rbind,iris.apply)
  ven2 = t(iris.combine) %>% as.data.frame()

  head(ven2)
  A <- combn(colnames(ven2),3)
  ven2$mean = rowMeans(ven2)
  tax = ggClusterNet::vegan_tax(ps)
  otutax = cbind(ven2,tax)
  head(otutax)
  otutax$Phylum[otutax$Phylum == ""] = "Unknown"


  # i = 1

  plot = list()
  for (i in 1:dim(A)[2]) {
    x = A[1,i]
    y = A[2,i]
    z = A[3,i]
    p <- ggtern::ggtern(data=otutax,aes_string(x = x,y=y,z=z,color = "Phylum",size ="mean" ))+
      geom_point() + theme_void()
    p
    plot[[i]] = p
    names(plot)[i] =  paste(paste(x,y,z,sep = "_"),sep = "")


  }
  return(list(plot,dataplot = otutax,groups = A))
}




