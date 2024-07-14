remove_rankID2 = function(ps) {
  tax = ps %>% vegan_tax() %>% as.data.frame()
  print(head(tax))
  tax$subkingdom =  NULL
  tax$Kingdom = tax$Kingdom %>% strsplit( "[__]") %>%
    sapply(`[`, 3)

  tax$Phylum = tax$Phylum %>% strsplit( "[__]") %>%
    sapply(`[`, 3)
  tax$Class = tax$Class %>% strsplit( "c__") %>%
    sapply(`[`, 2)
  tax$Class = tax$Class %>% strsplit( "_[a-z]__") %>%
    sapply(`[`, 1)

  tax$Order = tax$Order %>% strsplit( "o__") %>%
    sapply(`[`, 2)
  tax$Order  = tax$Order  %>% strsplit( "_[a-z]__") %>%
    sapply(`[`, 1)

  tax$Family = tax$Family %>% strsplit( "f__") %>%
    sapply(`[`, 2)
  tax$Family  = tax$Family  %>% strsplit( "_[a-z]__") %>%
    sapply(`[`, 1)


  tax$Genus = tax$Genus %>% strsplit( "g__") %>%
    sapply(`[`, 2)
  tax$Genus  = tax$Genus  %>% strsplit( "_[a-z]__") %>%
    sapply(`[`, 1)
  tax$Species = tax$Species %>% strsplit( "s__") %>%
    sapply(`[`, 2)
  tax$Species  = tax$Species  %>% strsplit( "_[a-z]__") %>%
    sapply(`[`, 1)
  print(head(tax))
  tax_table(ps) = tax_table(as.matrix(tax))
  return(ps)
}


change.OTU.name.ps2 = function(ps0 = ps,newnm ="metab_id" ){
  otu = ps0 %>% vegan_otu() %>% t() %>%
    as.data.frame()
  head(otu)
  rep = row.names(otu)
  row.names(otu) = tax[[newnm]]

  tax = ps0 %>% vegan_tax() %>%
    as.data.frame()
  tax$oldnm = row.names(tax)

  row.names(tax) = tax[[newnm]]
  head(tax)

  ps = phyloseq(
    otu_table(as.matrix(otu),taxa_are_rows =TRUE),
    tax_table(as.matrix(tax)),
    sample_data(ps0)
  )
  return(ps)
}


my.env = function(){
  library(EasyMultiOmics)
  library(ggClusterNet)
  library(phyloseq)
  library(tidyverse)
}


remove.all.zone = function(ps,n = 0){
  ps1 = ps %>%  filter_taxa(function(x) sum(x ) > n , TRUE)
  return(ps1)
}
