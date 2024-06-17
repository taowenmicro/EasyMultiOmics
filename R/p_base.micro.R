p_base.micro = function(ps,Top = 300) {
  alltax = ps %>%
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

