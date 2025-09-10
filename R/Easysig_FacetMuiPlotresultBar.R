#' Use facets to show the results of variance analysis, bar plot selected
#'
#' @param data a data.frame contain the input data;Enter the data frame, the first column is the sample number, the second column is the group, note that the group label must be set to group, and the third column and later are the indicators for measurement or collection.
#' @param  num col index wtich need to test
#' @param result output from aovMcomper or KwWlx. You can also import result calculated from other software (a data frame)
#' @param sig_show Distinctive display, "abc" or "line"
#' @param ncol If you choose faceted display, choose to place several graphics per line
#' @examples
#' # data(data_wt)
#' result = MuiKwWlx(data = data_wt,num = c(4:6))
#' result1 = FacetMuiPlotresultBar(data = data_wt,num = c(4:6),result = result,sig_show ="abc",ncol = 2 )
#' result1[[1]]
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export




FacetMuiPlotresultBar = function(data = data_wt,num = c(4:6),result = result,sig_show ="abc",ncol = 3,
                                 fac.level  = NULL,
                                 mult.y = 0.3
                                 ){
  N = num[1]
  name = colnames(data[N])

  as = result[match( name,colnames(result))]
  colnames(as) = "groups"
  as$group = row.names(as)

  PlotresultBar = aovMuiBarPlot(data = data, i= N,sig_show =sig_show,result = as)

  p = PlotresultBar[[2]]
  p
  name = colnames(data[N])
  p$name = name
  A = p

  for (N in num[-1]) {

    name = colnames(data[N])

    as = result[match( name,colnames(result))]
    colnames(as) = "groups"
    as$group = row.names(as)

    PlotresultBox = aovMuiBarPlot(data = data, i= N,sig_show =sig_show,result = as)

    p = PlotresultBox[[2]]

    name = colnames(data[N])
    p$name = name

    A = rbind(A,p)
  }
  # head(A)

  # a = max(A$SD)*1.2
  if (!is.null(fac.level)) {
  A$name  = factor(A$name,levels = fac.level)
  }
  lab = A %>% distinct(group ,name,.keep_all = TRUE)
  p<-ggplot(A, aes(x=group , y=mean ))+
    geom_bar(aes(fill = group),stat = "identity", width = 0.4,position = "dodge") +
    geom_bar(data = A, aes(x=1 , y= (mean + SD)*1.1),stat = "identity", width = 0.4,position = "dodge",alpha = 0) +
    geom_errorbar(aes(ymin=mean - SD,
                      ymax=mean+SD),
                  colour="black",width=0.1,linewidth = 0.5)+

    scale_y_continuous(expand = expansion(mult = c(0.05, mult.y)))+
    labs(
      x= NULL,
      y=NULL
      # title = paste("Normality test",p1,"Homogeneity of variance",p2,sep = ":")
    ) +

    geom_text(data=lab, aes(x=group, y=mean + SD,label=groups),vjust = -1,alpha = 0.6)+
    guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL)) +
    facet_wrap(.~name,scales="free_y",ncol  = ncol)


  p

  if (length(unique(data$group))>3){	p=p+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))}
  p


  return(list(p,table = A))
}



