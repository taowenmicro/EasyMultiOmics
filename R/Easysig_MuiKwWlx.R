#' Multi-group data for non-parametric test
#'
#' @param data a data.frame contain the input data
#' @param num index col which need plot
#' @examples
#' # data(data_wt)
#'  result = MuiKwWlx(data = data_wt,num = c(4:6))
#' result
#' @return data frame
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


#Multi-group data for non-parametric test
MuiKwWlx = function(data = data_wt,num = c(4:6)){
  N = num[1]

  data_wt = data
  result = KwWlx(data = data_wt, i= N)
  aa = result[[1]]
  name = colnames(data_wt[N])

  colnames(aa)[1] = name
  aa$group = NULL
  A = aa

  for (N in num[-1]) {
    result = KwWlx(data = data_wt, i= N)
    aa = result[[1]]
    name = colnames(data_wt[N])

    colnames(aa)[1] = name
    aa$group = NULL

    A <- cbind(A,aa)

  }

  return(A)

}
