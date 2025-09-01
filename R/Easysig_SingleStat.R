#' A complete single-factor difference analysis process for single col data, including normal test, homogeneity analysis of variance, selection of variance analysis or non-parametric test, selection of candidate d visualization method.
#'
#' @param data a data.frame contain the input data.Enter the data frame, the first column is the sample number,
#'  the second column is the group, note that the group label must be set to group, and the third column and later are the
#'  indicators for measurement or collection.
#' @param  num col index wtich need to test
#' @param method_cv which method to test  variance, included "leveneTest","bartlett.test"
#' @param method_Mc which specify the method of multiple comparisons.included:"Tukey","LSD","SNK","Duncan","scheffe","scheffe".
#' @param sig_show Distinctive display, "abc" or "line"
#' @examples
#' # data(data_wt)
#' result = SingleStat(data = data_wt,plot = "bar",method_Mc = "Tukey",i= 4,sig_show ="abc")
#' result[[1]]
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


SingleStat = function(data = data_wt,plot = "bar",method_Mc = "Tukey",i= 4,sig_show ="abc"){

  NorCV = NorNorCVTest(data = data, i= i,method_cv = "leveneTest")
  #
  a = NorCV[[1]]
  p1 = length(a$p.value[-dim(a)[1]][!a$p.value[-dim(a)[1]] >= 0.05]) == 0
  p2 = NorCV[[2]]


  if (p1!= 0& p2 >.05) {
    result= aovMcomper2 (data = data, i= i,method_Mc = "Tukey")
    A = print("aov")
    A
  } else if (p1  != 0| p2 <.05){
    result = KwWlx2(data = data, i= i)
    A = print("wlx")
    A
  }

  if (plot == "bar") {
    PlotresultBar = aovMuiBarPlot(data = data, i= i,sig_show =sig_show,result = result[[1]])
    p = PlotresultBar[[1]]

  } else if (plot == "box"){
    PlotresultBox = aovMuiBoxP2(data = data, i= i,sig_show =sig_show,result = result[[1]])
    p = PlotresultBox[[1]]

  }else if (plot == "boxbar"){
    PlotresultBox = aovMuiBoxBarP(data = data, i= i,sig_show =sig_show,result = result[[1]])
    p = PlotresultBox[[1]]

  }

  return(list(p,table = result[[1]],method =A))

}

