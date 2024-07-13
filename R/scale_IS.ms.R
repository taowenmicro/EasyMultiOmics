

# ps.ms2 = scale_IS.ms(ps = ps.ms,IS = "metab_8497")
scale_IS.ms = function(
    ps = ps.ms,
    IS = "metab_8497"
){
  otu = ps %>% vegan_otu() %>% t()
  head(otu)

  dat = otu[IS,]
  i = 1
  for (i in 1:dim(otu)[2]) {
    otu[,i] = otu[,i]/dat[i]
  }
  otu2 = otu1[row.names(otu1) != IS,]
  otu_table(ps) = otu_table(as.matrix(otu2),taxa_are_rows = TRUE)
  return(ps)
}
