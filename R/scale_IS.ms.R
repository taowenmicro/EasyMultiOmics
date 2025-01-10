#' @title Scale metabolome data by normalizing based on selected metabolite(s)
#' @description
#' The scale_IS.ms function scales  metabolite data by normalizing based on selected metabolite(s).
#' It  standardizes the  metabolite data of each metabolite  by the expression levels of selected metabolite(s) in each sample.
#' @param ps A phyloseq format file used as an alternative for the input containing metabolite, tax, and sample data.
#' @param IS A character vector specifying the selected metabolite(s).
#' @return The scaled phyloseq format file after normalization based on the selected metabolite(s).
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' ps.ms2 = scale_IS.ms(ps = ps.ms,IS = "metab_8497")
#' @export

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
  otu2 = otu[row.names(otu) != IS,]
  otu_table(ps) = otu_table(as.matrix(otu2),taxa_are_rows = TRUE)
  return(ps)
}
