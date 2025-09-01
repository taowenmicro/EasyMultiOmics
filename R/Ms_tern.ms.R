#' @title Generate Ternary Plots for metabolites composition data
#'
#' @description
#' This function generates ternary plots for visualizing metabolites composition data in three groups.
#' It uses relative abundance data from a `phyloseq` object and displays metabolites distribution across specified groupings.
#'
#' @param ps A `phyloseq` object containing metabolites composition data.
#'
#' @return
#' A list containing the following elements:
#' \describe{
#'   \item{plot}{A list of ternary plots (`ggtern` objects) for different group combinations.}
#'   \item{dataplot}{A data frame with combined metabolites abundance and taxonomic information.}
#'   \item{groups}{A matrix showing all possible combinations of three groups used for ternary plots.}
#' }
#'
#' @details
#' This function performs the following steps:
#' \itemize{
#'   \item Standardizes the metabolites composition table to relative abundances using `phyloseq::transform_sample_counts`.
#'   \item Splits the OTU data into groups based on the metadata column `Group`.
#'   \item Computes the mean abundance for each metabolite within each group.
#'   \item Generates ternary plots for all possible combinations of three groups using the `ggtern` package.
#' }
#'
#' The plots visualize the relative abundances of metabolites in three groups, with points sized by their mean abundance and colored by their taxonomic classification (mode).
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' \dontrun{
#' ps1 = ps.ms %>% filter_OTU_ps(500)
#' res = Micro_tern.ms(ps1)
#' p15 = res[[1]]
#' dat =  res[[2]]
#' }
#'
#' @export
Ms_tern.ms = function (ps = ps,color = color)
{
  ps_rela = phyloseq::transform_sample_counts(ps, function(x) x/sum(x))
  ps_rela
  otu = ggClusterNet::vegan_otu(ps_rela) %>% as.data.frame()
  iris.split <- split(otu, as.factor(as.factor(phyloseq::sample_data(ps)$Group)))
  iris.apply <- lapply(iris.split, function(x) colSums(x[]))
  iris.combine <- do.call(rbind, iris.apply)
  ven2 = t(iris.combine) %>% as.data.frame()
  head(ven2)
  A <- combn(colnames(ven2), 3)
  ven2$mean = rowMeans(ven2)
  tax = ggClusterNet::vegan_tax(ps)
  otutax = cbind(ven2, tax)
  head(otutax)

  otutax
  plot = list()
  for (i in 1:dim(A)[2]) {
    #i=1
    x = A[1, i]
    y = A[2, i]
    z = A[3, i]
    colnames(otutax)
    p <- ggtern::ggtern(data = otutax, aes_string(x = x,
                                                  y = y,
                                                  z = z, color = color,
                                                  size = "mean")) +
      geom_point() + theme(legend.position = "none")
    p
    plot[[i]] = p
    names(plot)[i] = paste(paste(x, y, z, sep = "_"), sep = "")
  }
  return(list(plot, dataplot = otutax, groups = A))
}
