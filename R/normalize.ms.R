#' @title Normalize metabolite data using specified method
#' @description
#' The normalize.ms function normalizes  metabolite data based on the specified method.
#' Three normalization methods are supported: "rela" for relative abundance normalization,
#' "sampling" for sampling normalization (adjusting counts to make the total count of each sample equal to the average total count),
#' and "log" for log transformation.
#' @param ps A phyloseq format file used as an alternative for the input containing metabolite, tax, and sample data.
#' @param method The normalization method to be applied,default is "rela".
#' @return The new phyloseq format file by Selected method after the standardization.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' ps.ms2 = normalize.ms(ps = ps.ms,method = "rela")
#' @export

normalize.ms = function(ps,
                        method = "rela"
){
  if (method == "rela") {
    ps1  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x,na.rm =TRUE) )

  }

  if (method == "sampling" ) {
    total = mean(sample_sums(ps));total
    standf = function(x,t = total)round(t*(x/sum(x)))
    ps1 = phyloseq::transform_sample_counts(ps,standf)
  }
  if (method == "log") {
    ps1  = phyloseq::transform_sample_counts(ps, function(x) log(1 + x) )

  }
  return(ps1)
}
