#' @title Scale  metabolite data by normalizing based on Quality Control (QC) samples
#' @description
#' The scale_QC.ms function scales metabolite data by normalizing based on Quality Control (QC) samples.
#' It calculates the average value of each metabolite in the QC samples and standardizes  metabolite data of each sample by this average value.
#' @param ps A phyloseq format file used as an alternative for the input containing metabolite, tax, and sample metadata.
#' @param QC A character vector specifying the QC sample names.
#' @return The scaled phyloseq format file after normalization based on the QC samples.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' ps.ms2 = scale_QC.ms(ps = ps.ms,QC = c("OE36_2","OE36_3","OE2_1"))
#' @export


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


