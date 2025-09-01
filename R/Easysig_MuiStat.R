#' A complete single-factor difference analysis process for mui col data, including normal test, homogeneity analysis of variance, selection of variance analysis or non-parametric test, selection of candidate d visualization method.
#'
#' @param data a data.frame contain the input data;Enter the data frame, the first column is the sample number, the second column is the group, note that the group label must be set to group, and the third column and later are the indicators for measurement or collection.
#' @param  num col index wtich need to test
#' @param method_cv which method to test  variance, included "leveneTest","bartlett.test"
#' @param method_Mc which specify the method of multiple comparisons.included:"Tukey","LSD","SNK","Duncan","scheffe","scheffe".
#' @param sig_show Distinctive display, "abc" or "line"
#' @param plot  chose the plot type shuld be a barplot or boxplot or bar+box plot
#' @param plottype select a plot type, included "single" and "mui"
#' @param ncol If you choose faceted display, choose to place several graphics per line
#' @param path creat a folder to save the single plot, if you seleted the model "single"
#' @examples
#' # data(data_wt)
#' result = MuiStat(data = data_wt,num = c(4,5,6),method_cv = "leveneTest",method_Mc = "Tukey",sig_show  = "abc",ncol = 2,plot = "boxbar",plottype = "mui")
#' result[[1]]
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export



MuiStat = function(data = data_wt,num = c(4,5,6),method_cv = "leveneTest",method_Mc = "Tukey",sig_show  = "abc",ncol = 2,plot = "bar",plottype = "mui",path = "./output"){

  norCv = MuiNorCV(data = data,num = num,method_cv = "leveneTest")

  # Divide the group that conforms to the normal distribution and the homogeneity of the variance into one group,
  # and divide the group that does not meet the two into a group;
  AA = c()
  BB = c()
  for (i in 1:length(num)) {

    if (norCv[,"cor"][i] == TRUE & norCv[,"CV"][i] == TRUE) {
      num[i]
      AA = c(AA,num[i])
      AA

    }

    if (norCv[,"cor"][i] == FALSE | norCv[,"CV"][i] == FALSE) {
      num[i]
      BB = c(BB,num[i])
      BB

    }

  }

  if (!is.null(AA)) {
    resultAA = MuiaovMcomper2(data = data,num = AA,method_Mc = method_Mc )
  } else if(is.null(AA)){
    resultAA = NULL

  }

  if(!is.null(BB)){
  resultBB = MuiKwWlx2(data = data,num = BB)
  } else if(is.null(BB)){
    resultBB = NULL

  }

  if (!is.null(AA) & !is.null(BB)) {
    resultall = cbind(resultAA,resultBB)
  } else if (is.null(AA)) {
    resultall = cbind(resultBB)
  } else if (is.null(BB)) {
    resultall = cbind(resultAA)
  }




  num = c(AA,BB)
  if (plottype == "single") {
    if (plot == "bar") {
      plot = MuiPlotresultBar(data = data,num = num,result = resultall,sig_show = sig_show,path = path)
      p = "Folder"
    }
    if (plot == "box") {
      plot = MuiPlotresultBox(data = data,num = num,result = resultall,sig_show = sig_show,path = path)
      p = "Folder"
    }
    if (plot == "boxbar") {
      plot = MuiPlotresultBox(data = data,num = num,result = resultall,sig_show = sig_show,path = path)
      p = "Folder"
    }
  }
  if (plottype == "mui") {
    if (plot == "bar") {
      result1 = FacetMuiPlotresultBar(data = data,num = num,result = resultall,sig_show =sig_show,ncol = ncol )
      p = result1[[1]]
    }
    if (plot == "box") {
      result1 = FacetMuiPlotresultBox(data = data,num = num,result = resultall,sig_show =sig_show,ncol =ncol )
      p = result1[[1]]
    }
    if (plot == "boxbar") {
      result1 = FacetMuiPlotReBoxBar(data = data,num = num,result = resultall,sig_show =sig_show,ncol =ncol )
      p = result1[[1]]
    }
  }

  return(list(p, aov = AA,wlx = BB,table =resultall ))


}
