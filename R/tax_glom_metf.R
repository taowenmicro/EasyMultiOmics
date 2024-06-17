tax_glom_metf <- function(ps = ps,ranks = "level1") {


  if (  is.numeric(ranks)) {
    ranks <- phyloseq::rank.names(ps)[ranks]
  }


  otu <- as.data.frame(t(vegan_otu(ps)))
  tax <- as.data.frame(vegan_tax(ps))
  head(tax)
  # building group
  tax[[ranks]][is.na(tax[[ranks]])] = "Unknown"
  tax[[ranks]][tax[[ranks]] == ""] = "Unknown"
  tax[[ranks]][tax[[ranks]] == "NA"] = "Unknown"
  split <- split(otu,tax[[ranks]])
  #calculate sum by group
  apply <- lapply(split,function(x)colSums(x[]))
  # result chack
  otucon <- do.call(rbind,apply)

  head(tax)
  # row.names(otucon) %>% unique()
  # tax$level2 %>% unique()
  # tax[[ranks]] %>% unique()

  taxcon <- tax[1:match(ranks,colnames(tax))] %>% as.matrix()
  taxcon <- taxcon[!duplicated(tax[[ranks]]),] %>% as.data.frame()
  head(taxcon)
  # colnames(taxcon) = j
  if (ncol(taxcon) == 1) {
    colnames(taxcon) = j
  }
  #-tax name with NA wound be repeated with unknown
  taxcon[[ranks]][is.na(taxcon[[ranks]])] = "unknown"
  row.names(taxcon) <- taxcon[[ranks]]

  match(row.names(otucon),row.names(taxcon))
  head(tax)
  pscon <- phyloseq::phyloseq(
    phyloseq::otu_table( as.matrix(otucon),taxa_are_rows = TRUE),
    phyloseq::tax_table(as.matrix(taxcon)),
    phyloseq::sample_data(ps)
  )
  return(pscon)
}
