# ps.ms2 = scale_QC.ms(
#   ps = ps.ms,
#   QC = c("OE36_2","OE36_3","OE2_1"))


scale_QC.ms = function(
    ps = ps.ms,
    QC = c("OE36_2","OE36_3","OE2_1")
){
  otu = ps %>% vegan_otu() %>% t() %>% as.data.frame()
  tem = otu[,QC]
  tem2 = rowSums(tem)
  for (i in 1:length(tem2)) {
    otu[i,] = otu[i,]/tem2[i]
  }

  otu2 = otu[,!colnames(otu) %in% QC]
  otu_table(ps) = otu_table(as.matrix(otu2),taxa_are_rows = TRUE)
  return(ps)
}


