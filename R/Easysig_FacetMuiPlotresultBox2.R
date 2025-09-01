#' Use facets to show the results of variance analysis, box plot selected
#'
#' @param data a data.frame contain the input data;Enter the data frame, the first column is the sample number, the second column is the group, note that the group label must be set to group, and the third column and later are the indicators for measurement or collection.
#' @param  num col index wtich need to test
#' @param result output from aovMcomper or KwWlx. You can also import result calculated from other software (a data frame)
#' @param sig_show Distinctive display, "abc" or "line"
#' @param ncol If you choose faceted display, choose to place several graphics per line
#' @param loc cahnge the label of "abc" to Suitable location
#' @examples
#' # data(data_wt)
#' result = MuiKwWlx(data = data_wt,num = c(4:6))
#' result1 = FacetMuiPlotresultBox(data = data_wt,num = c(4:6),result = result,sig_show ="abc",ncol = 2 )
#' result1[[1]]
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{2018203048@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export

FacetMuiPlotresultBox2 = function(data = data_wt,num = c(4:6),result = result,sig_show ="abc",ncol = 3,fac.level  = NULL ){
  for (N in num) {

    name = colnames(data[N])

    as = result[match( name,colnames(result))]
    # as = result[c(N - length(num))]
    as
    colnames(as) = "groups"
    as$group = row.names(as)

    PlotresultBox = aovMuiBoxP2(data = data, i= N,sig_show =sig_show,result =as)
    PlotresultBox[[1]]
    p = PlotresultBox[[2]]
    p
    name = colnames(data[N])
    p$name = name

    if (N == num[1]) {
      A = p
    }
    if (N != num[1]) {
      A = rbind(A,p)
    }

  }

  if (!is.null(fac.level)) {
    A$name  = factor(A$name,levels = fac.level)
  }
  p<-ggplot(A, aes(x=group , y=dd ))+ geom_boxplot(alpha=1, aes(fill=group)) +
    geom_jitter( position=position_jitter(0.17), size=0.1, alpha=0.5)+
    labs(x="", y="")+
    # theme_classic()+
    geom_text(data=A, aes(x=group , y=y ,label=stat))+
    #scale_colour_manual(values = mi)+
    #scale_fill_manual(values = mi)+
    #labs(title = "toamto hea and dis")+
    guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+facet_wrap(.~name,scales="free_y",ncol  = ncol)
  p

  p <- ggpubr::ggpaired(A, x = "group", y = "dd",id="pid",
                  color = "group", palette = "lancet", point.size=2,linetype=2,
                  line.color = "gray", outlier.size=0, size=0.7, width=0.5, fill="transparent",
                  short.panel.labs = FALSE) +
    labs(x="", y="")+
    # theme_classic()+
    geom_text(data=A, aes(x=group , y=y ,label=stat))+
    #scale_colour_manual(values = mi)+
    #scale_fill_manual(values = mi)+
    #labs(title = "toamto hea and dis")+

    facet_wrap(.~name,scales="free_y",ncol  = ncol)





  return(list(p,table = A))
}
