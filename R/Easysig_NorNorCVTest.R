#' Test normality and homogeneity of variance of Input data
#'
#' @param data a data.frame contain the input data
#' @param  i col index wtich need to test
#' @param method_cv which method to test  variance, included "leveneTest","bartlett.test".
#' @examples
#' # data(data_wt)
#' NorNorCVTest(data = data_wt, i= 4)
#' @return data frame
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


NorNorCVTest = function(data = data_wt, i= 4,method_cv = "leveneTest",...){
  ##----Test normality and homogeneity of variance of Input data-------
  ss <- data %>%
    dplyr::select("group",count = i)
  shapiro.test.multi <- function(
    data,
    value,
    group)
  {
    table(data[,group]) %>%
      data.frame(.) -> a1
    a2 <- as.vector(a1[,1])
    data = data.frame(group = data[,group],  #
                      value = data[,value])  #

    test.result <- data.frame(No=0,        #row number
                              Name=0,      #group name
                              W=0,         #W value
                              p.value=0,   #p value
                              norm.test=0) #result
    for (i in (1:length(a2))){
      # subgroup for shapiro.text
      subset(data,
             group == a2[i],
             select = value) %>%
        .[,1] %>%
        shapiro.test(.) -> t.r
      test.result[i,1] = i              #group number
      test.result[i,2] = a2[i]          #group name
      test.result[i,3] = t.r$statistic  #w value
      test.result[i,4] = t.r$p.value    #p value
      if
      (t.r$p.value > 0.05)
        test.result[i,5] = "Norm"
      else
        test.result[i,5] = "Other_situation"
    }
    test.result[nrow( test.result)+1,1] = "Test Method:"  # final col add  the method of text
    test.result[nrow( test.result),2] = "Shapiro-Wilk"    #final col two row add the method
    return(test.result)
  }
  a = shapiro.test.multi(data[c(2,i)],value = colnames(data[c(2,i)][2]),group = "group")
  a
  # selected the method for tast the Homogeneity of variance
  # p2 >=.05:Homogeneity of variance
  if (method_cv == "leveneTest" ) {
    xc <- car::leveneTest(count~group,data=ss)
    p2 <- xc[[3]][1]
    p2 <- round(p2,3)
  }
  if (method_cv == "bartlett.test" ) {
    xc <- bartlett.test(count~group,data=ss)
    p2 <- xc[[3]]
    p2 <- round(p2,3)
  }
  return(list(a,p2))
}
