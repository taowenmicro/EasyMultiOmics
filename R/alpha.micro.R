

#' Title Calculate the diversity for each sample or group
#'
#' @param otu OTU/ASV table;
#' @param tax taxonomy table;
#' @param map matrix or data frame, including ID and group;
#' @param ps phyloseq class;
#' @param group column name for group, such as "Group".
#' @param sampling sampling OTU/ASV table with the minisize sequence count;
#'
#' @return data frame of sample diversity indices.
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @export
#'
#' @examples
alpha.micro = function(otu = NULL,
                 tax = NULL,
                 map = NULL,
                 ps = NULL,
                 group = "Group",

                 sampling = TRUE){

  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)
  if (sampling == TRUE) {
    samplesize = min(phyloseq::sample_sums(ps))
    if (samplesize == 0) {
      print("0 number sequence of some samples")
      print("median number were used")
      ps11  = phyloseq::rarefy_even_depth(ps,sample.size = samplesize)
    } else{
      ps11  = phyloseq::rarefy_even_depth(ps,sample.size = samplesize)
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






