#' @title Visualize taxonomic information using circular tree layout
#' @description
#' convert dataframe contained hierarchical relationship or other classes to
#' treedata class and visualize the circular tree.
#' @param ps A phyloseq format file used as an alternative for the input containing otu, tax, and map.
#' @param Top The top microorganisms to convert and visualize.
#' @returns A ggplot object displaying taxonomic information in a circular tree layout.
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' p1 <- p_base.micro(ps.16s,Top = 100)
#' p1
p_base.micro = function(ps,Top = 300) {
  alltax = ps.16s %>%
    ggClusterNet::filter_OTU_ps(Top) %>%
    ggClusterNet::vegan_tax() %>%
    as.data.frame()
  alltax$OTU = row.names(alltax)
  alltax$Kingdom = paste(alltax$Kingdom,sep = "_Rank_")
  alltax$Phylum = paste(alltax$Kingdom,alltax$Phylum,sep = "_Rank_")
  alltax$Class = paste(alltax$Phylum,alltax$Class,sep = "_Rank_")
  alltax$Order = paste(alltax$Class,alltax$Order,sep = "_Rank_")
  alltax$Family = paste(alltax$Order,alltax$Family,sep = "_Rank_")
  alltax$Genus = paste(alltax$Family,alltax$Genus,sep = "_Rank_")
  alltax$Species = paste(alltax$Genus,alltax$Species,sep = "_Rank_")
  alltax[is.na(alltax)] = "Unknown"
  trda <- MicrobiotaProcess::convert_to_treedata(alltax)

  p <- ggtree(trda, layout="circular", size=0.2, xlim=c(30,NA)) +
    geom_point(
      pch = 21,
      size=3,
      alpha=1,
      fill = "#FFFFB3"
    )
  p$data$lab2 <- p$data$label %>% strsplit( "_Rank_") %>%
    sapply(
      function(x) x[length(x)]
    )
  p$data$lab2   = gsub("st__","",p$data$lab2 )
  p$data$nodeSize = 1
  return(p)
}

