
#' @title Calculate Phylogenetic Diversity (PD) for Samples
#' @description
#'This function calculates the phylogenetic diversity (PD) for each sample in a given phyloseq object.
#' @param ps A `phyloseq` object. It must contain an OTU/ASV table (`otu_table`), a phylogenetic tree (`phy_tree`), and sample metadata (`sample_data`).
#' @param group column name for group, such as "Group".
#'
#' @return data frame of sample pd indices.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @export
#'
#' @examples
#' \dontrun{
#' library(ape)
#' library(picante)
#' tab2 = alpha.pd(ps.16s)
#' head(tab2)
#' }
alpha.pd = function(ps=NULL,group = "Group"){
  com_2020 <- ps %>% vegan_otu() %>%
    as.data.frame()

  rooted <- phy_tree(ps)

  cover2020.pd<-pd(com_2020,rooted,include.root=F)

  map = sample_data(ps)
  head(map)

  data = cbind(map[,c("ID","Group")],pd = cover2020.pd[,1])
  head(data)
  colnames(data)[2] = "group"
  data$group = as.factor(data$group)
  return(data)
}
