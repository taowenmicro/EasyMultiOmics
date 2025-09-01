#' Perform ANOVA and run multiple comparisons and order in order to mean value of treatment
#'
#' @param data a data.frame contain the input data
#' @param  i col index wtich need to test
#' @param method_Mc which specify the method of multiple comparisons.included:"Tukey","LSD","SNK","Duncan","scheffe","scheffe".
#' @examples
#' # data(data_wt)
#' result= aovMcomper2 (data = data_wt, i= 5,method_Mc = "Tukey")
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
aovMcomper2 = function( data = data_wt, i= 3,method_Mc = "Tukey"){
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
    aa <- insx$mcletters$monospacedLetters %>% as.data.frame()
    colnames(aa) = c("groups")
    aa$groups = gsub(" ","",aa$groups)
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

  dat <- ss %>% group_by(group) %>%
    dplyr::summarise(mean = mean(count)) %>%
    dplyr::inner_join(as_tibble(aa),by = c("group" = "group")) %>%
    dplyr::arrange(desc(mean))
  tmp.1 <- dat$groups %>% unique()
  tmp.2 <- data.frame(ori = tmp.1,new = sort(tmp.1,decreasing = FALSE))
  i = 1
  for (i in 1:nrow(dat)) {
    dat[i,3] <- tmp.2$new[match(dat[i,3],tmp.2$ori)]
  }
  dat$mean = NULL
  dat = dat[match(dat$group,aa$group),][,c(2,1)] %>% as.data.frame()
  row.names(dat) = dat$group
  # aa =ord_sig(data = aa,ID = "groups")
  return(list(Muicomper = dat,model))
}

#'  Easy for Non-parametric test and order in order to mean value of treatment
#'
#' @param data a data.frame contain the input data
#' @param  i col index wtich need to test
#' @param method wilcox.test or t.test could be seleced.
#' @examples
#' # data(data_wt)
#' result = KwWlx2(data = data_wt, i= 4)
#' result[[1]]
#' @return list with two data frame
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{2018203048@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export

KwWlx2 = function(data = data_wt, i= 3,method = "wilcox.test"){

  ss <- data %>%
    dplyr::select("group",count = i)

  # kruskal.test
  krusk = ggpubr::compare_means(count ~ group, data=ss, method = "kruskal.test")
  # krusk = kruskal.test(count ~ group, data=ss)

  sumkrusk=as.data.frame(krusk)
  sumkrusk
  #<0.05,It shows that there are differences between multiple groups, you can conduct pairwise non-parametric tests, and mark the letters
  krusk <- ggpubr::compare_means(count ~ group, data=ss, method = method)
  xx=as.data.frame(krusk)

  #mean for order
  da <- ss %>%
    dplyr::group_by(group) %>%
    dplyr::summarise( mean = mean(count))  %>%
    dplyr::arrange(desc(mean)) %>%
    as.data.frame()

  if (as.character(da$group)[1]  %in% xx$group1) {
    xx = xx[order(factor(xx$group1,as.character(da$group))),]
    wilcox_levels = xx$p
    names(wilcox_levels) =  paste(xx$group1,xx$group2,sep = "-")
  } else {
    xx = xx[order(factor(xx$group2,as.character(da$group))),]
    wilcox_levels = xx$p
    names(wilcox_levels) =  paste(xx$group2,xx$group1,sep = "-")
  }


  wilcox.labels <- data.frame(multcompView::multcompLetters(wilcox_levels, threshold = 0.05)['Letters'])
  colnames(wilcox.labels) = "groups"
  aa = wilcox.labels
  aa$group = row.names(aa)
  # aa =ord_sig(data = aa,ID = "groups")
  dat <- ss %>% group_by(group) %>%
    dplyr::summarise(mean = mean(count)) %>%
    dplyr::inner_join(as_tibble(aa),by = c("group" = "group")) %>%
    dplyr::arrange(desc(mean))
  tmp.1 <- dat$groups %>% unique()
  tmp.2 <- data.frame(ori = tmp.1,new = sort(tmp.1,decreasing = FALSE))
  i = 1
  for (i in 1:nrow(dat)) {
    dat[i,3] <- tmp.2$new[match(dat[i,3],tmp.2$ori)]
  }
  dat$mean = NULL
  dat = dat[match(dat$group,aa$group),][,c(2,1)] %>% as.data.frame()
  row.names(dat) = dat$group
  return(list(dat,wilcox = krusk,kruskal = sumkrusk))
}


#' Perform ANOVA and run multiple comparisons for more sets of data and order in order to mean value of treatment
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

MuiaovMcomper2 = function(data = data_wt,num = c(4:6),method_Mc = "Tukey"){
  data$group = as.factor(data$group)
  N = num[1]
  result = aovMcomper2 (data = data, i= N,method_Mc = method_Mc)
  aa = result[[1]]
  name = colnames(data[N])
  colnames(aa)[1] = name
  aa$group = NULL
  A = aa

  for (N in num[-1]) {
    result = aovMcomper2 (data = data, i= N,method_Mc = method_Mc)
    aa = result[[1]]
    name = colnames(data[N])

    colnames(aa)[1] = name
    aa <- aa[match(row.names(A),row.names(aa)),]
    aa$group = NULL
    A <- merge(A,aa,by = "row.names",all = T)
    row.names(A) = A$Row.names
    A$Row.names = NULL
  }

  return(A)
}


#' Multi-group data for non-parametric test
#'
#' @param data a data.frame contain the input data and order in order to mean value of treatment
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

MuiKwWlx2 = function(data = data_wt,num = c(4:6)){
  N = num[1]

  data_wt = data
  result = KwWlx2(data = data_wt, i= N)
  aa = result[[1]]
  name = colnames(data_wt[N])

  colnames(aa)[1] = name
  aa$group = NULL
  A = aa

  for (N in num[-1]) {
    result = KwWlx2(data = data_wt, i= N)
    aa = result[[1]]
    name = colnames(data_wt[N])

    colnames(aa)[1] = name
    aa$group = NULL

    A <- merge(A,aa,by = "row.names",all = T)
    row.names(A) = A$Row.names
    A$Row.names = NULL



  }

  return(A)

}












