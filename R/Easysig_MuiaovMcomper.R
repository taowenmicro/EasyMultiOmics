#' Perform ANOVA and run multiple comparisons for more sets of data
#'
#' @param data a data.frame contain the input data
#' @param num index col which to test normality and homogeneity
#' @param method_Mc which specify the method of multiple comparisons.included:"Tukey","LSD","SNK","Duncan","scheffe","scheffe".
#' @examples
#' # data(data_wt)
#' result = MuiaovMcomper(data = data_wt,num = c(4:6),method_Mc = "Tukey")
#' # utput result for multiple comparison results
#' result
#' @return data frame
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{2018203048@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export



MuiaovMcomper = function(data = data_wt,num = c(4:6),method_Mc = "Tukey"){
  N = num[1]

  result = aovMcomper (data = data, i= N,method_Mc = method_Mc)
  aa = result[[1]]
  name = colnames(data[N])
  colnames(aa)[1] = name
  aa$group = NULL
  A = aa

  for (N in num[-1]) {
    result = aovMcomper (data = data, i= N,method_Mc = method_Mc)
    aa = result[[1]]
    name = colnames(data[N])

    colnames(aa)[1] = name
    aa <- aa[match(row.names(A),row.names(aa)),]
    aa$group = NULL
    A =  cbind(A,aa)
  }

  return(A)
}
