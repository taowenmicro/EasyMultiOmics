#' @title Calculate the diversity for each sample or group
#' @description
#' The alpha.metf function calculates alpha diversity indices for metagenome functional composition.
#' This includes measures such as Shannon diversity, Pielou's evenness, and richness estimates (e.g., Chao1, ACE).
#'
#' @param otu Metagenome functional composition table;
#' @param tax Metagenomic Functional taxonomy table;
#' @param map A data frame or matrix containing metadata for samples. Must include at least `ID` and grouping column(s).
#' @param ps A `phyloseq` object. It must contain a metagenome functional composition table (`otu_table`), and sample metadata (`sample_data`).
#' @param group Column name for group, such as "Group".
#' @param sampling Sampling metagenome functional composition table with the minisize sequence count;
#'
#' @return  A data frame containing alpha diversity indices for each sample, along with the sample metadata. The output includes the following indices:
#' \itemize{
#'   \item `Shannon`: Shannon diversity index.
#'   \item `Inv_Simpson`: Inverse Simpson diversity index.
#'   \item `Pielou_evenness`: Pielou's evenness index.
#'   \item `Simpson_evenness`: Simpson's evenness index.
#'   \item `Richness`: Observed richness (number of genes).
#'   \item `Chao1`: Chao1 richness estimator.
#'   \item `ACE`: ACE richness estimator.
#' }
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @export
#'
#' @examples
#' \dontrun{
#' ps =ps.kegg %>% filter_OTU_ps(Top = 1000)
#' tab = alpha.metf(ps = ps,group = "Group")
#' }
alpha.metf = function(otu = NULL,
                       tax = NULL,
                       map = NULL,
                       ps = NULL,
                       group = "Group",

                       sampling = TRUE,
                       Plot = TRUE){

  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)
  if (sampling == TRUE) {
    samplesize = min(phyloseq::sample_sums(ps))
    if (samplesize == 0) {
      print("0 number sequence of some samples")
      print("median number were used")
      # ps11  = phyloseq::rarefy_even_depth(ps,sample.size = samplesize)
      ps11 <- phyRare_auto(ps = ps, N = samplesize, drop_insufficient = TRUE, rngseed = 123)
    } else{
      # ps11  = phyloseq::rarefy_even_depth(ps,sample.size = samplesize)
      ps11 <- phyRare_auto(ps = ps, N = samplesize, drop_insufficient = TRUE, rngseed = 123)
    }
  } else if(sampling == FALSE){
    ps11 = ps
  }
  mapping = phyloseq::sample_data(ps11)
  ps11 = phyloseq::filter_taxa(ps11, function(x) sum(x ) >0 , TRUE); ps11
  head(mapping)
  colnames(mapping) = gsub(group,"AA", colnames(mapping))

  mapping$Group = mapping$AA
  mapping$Group = as.factor(mapping$Group)
  mapping$Group

  count = as.data.frame(t(ggClusterNet::vegan_otu(ps11)))

  alpha=vegan::diversity(count, "shannon")

  x = t(count)
  head(x)

  Shannon = vegan::diversity(x)
  Shannon
  Inv_Simpson <- vegan::diversity(x, index = "invsimpson")
  Inv_Simpson

  S <- vegan::specnumber(x);S
  S2 = rowSums(x>0)


  Pielou_evenness <- Shannon/log(S)
  Simpson_evenness <- Inv_Simpson/S

  est <- vegan::estimateR(x)
  est <- vegan::estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]

  report = cbind(Shannon, Inv_Simpson, Pielou_evenness, Simpson_evenness,
                 Richness, Chao1,ACE)
  head(report)

  index = merge(mapping,report , by="row.names",all=F)

  return(index)
}






