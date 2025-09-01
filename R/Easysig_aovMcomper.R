#' Perform ANOVA and run multiple comparisons
#'
#' @param data a data.frame contain the input data
#' @param  i col index wtich need to test
#' @param method_Mc which specify the method of multiple comparisons.included:"Tukey","LSD","SNK","Duncan","scheffe","scheffe".
#' @examples
#' # data(data_wt)
#' result= aovMcomper (data = data_wt, i= 5,method_Mc = "Tukey")
#' # utput result for multiple comparison results
#' result[[1]]
#' @return data frame
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


# if (is.numeric(i)) {
#   ss <- data %>%
#     dplyr::select("group",count = i)
#
# } else if (is.character(i)) {
#   ss <- data %>%
#     dplyr::select("group",count = i)
# }

aovMcomper = function( data = data_wt, i= 3,method_Mc = "Tukey"){
  ss <- data %>%
    dplyr::select("group",count = i)

  # variance analysis
  model<-aov(count ~ group, data= ss)
  wtx1 = summary(model)
  wtx2 = wtx1[[1]]
  wtx3 = wtx2[5]#
  # Tukey
  if (method_Mc == "Tukey") {

    litter.mc <- multcomp::glht(model, linfct = multcomp::mcp(group = 'Tukey'))
    insx = multcomp::cld(litter.mc)
    aa <- insx$mcletters$monospacedLetters
    aa = as.data.frame(aa)
    colnames(aa) = c("groups")
    head(aa)
    aa$group = row.names(aa)

  }

  # LSD
  if (method_Mc == "LSD") {

    out <- agricolae::LSD.test(model,"group", p.adj="none")
    aa = out$group#
    aa$group = row.names(aa)
    aa = aa[2:3]
    aa
  }
  #SNK method (Student-Newman-Keuls)The results are similar to LSD.test.
  if (method_Mc == "SNK") {
    out <- agricolae::SNK.test(model,"group")
    aa = out$groups# View label of each group

    aa$group = row.names(aa)
    stat = aa
    aa = aa[2:3]
  }

  #Duncan
  if (method_Mc == "Duncan") {

    out <-agricolae::duncan.test (model,"group")
    aa = out$groups#  View label of each group

    aa$group = row.names(aa)
    stat = aa
    aa
    aa = aa[2:3]
  }

  #Scheffe features：The number of samples in each group is equal or unequal, but it is more used if the number of samples in each group is not equal;
  if (method_Mc == "scheffe") {
    out <-agricolae::scheffe.test (model,"group")
    aa = out$groups#  View label of each group
    aa$group = row.names(aa)
    stat = aa
    aa = aa[2:3]
  }
  # aa =ord_sig(data = aa,ID = "groups")
  return(list(Muicomper = aa,model))
}



