#' Test normality and homogeneity of variance of inputmore sets of data
#'
#' @param data a data.frame contain the input data
#' @param num index col which to test normality and homogeneity
#' @param method_cv which method to test  variance, included "leveneTest","bartlett.test".
#' @examples
#' # data(data_wt)
#' norCv = MuiNorCV(data = data_wt,num = c(4,5,6),method_cv = "leveneTest")
#' norCv
#' @return data frame
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export

# # 使用案例
# norCv = MuiNorCV(data = data_wt,num = c(4,5,6),method_cv = "leveneTest")
# norCv


## ---------------------第一部分----------------多组数据正态分布和方差齐性分析
# num = c(4,6)
# N =4
MuiNorCV = function(data = data_wt,num = c(4:6),method_cv = "leveneTest"){
  data_wt = data
  s1 = rep("A",length(num))
  s2 = rep("A",length(num))
  s0 = rep("A",length(num))
  i = 1
  for (N in num) {

    resul = NorNorCVTest(data = data_wt, i= N ,method_cv = "leveneTest")
    a = resul[[1]]
    b = resul[[2]]
    name = colnames(data_wt[N])
    a1 = length(a$p.value[-dim(a)[1]][!a$p.value[-dim(a)[1]] >= 0.05]) == 0
    b1 = b >=.05
    s1[i] = a1
    s2[i] = b1
    s0[i] = name
    i = i+1
  }

  result = cbind(s0,s1,s2)
  colnames(result) = c("DI","cor","CV")
  return(result)

}
