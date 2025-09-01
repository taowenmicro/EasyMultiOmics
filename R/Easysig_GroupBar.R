#' Variance test or non-parametric test results visualization, using boxbar
#'
#' @param data a data.frame contain the input data
#' @param num index col which need plot
#' @param sig_show Distinctive display, "abc" or we could find in the future
#' @param result output from aovMcomper or KwWlx. You can also import result calculated from other software (a data frame)
#' @param ylimit Adjust the y-axis coordinate range
#' @examples
#' # data(data_wt)
#' result = MuiKwWlx(data = data_wt,num = c(4:6))
#' result = GroupBar(data = data_wt,num = c(4:6),result = result,sig_show ="abc")
#' # output result
#' result[[1]]
#'
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export

GroupBar = function(data = data_wt,num = c(4:6),result = result,sig_show ="abc",ylimit = 1.2 ,
                    cord = FALSE
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
    p
    name = colnames(data[N])
    p$name = name

    A = rbind(A,p)
  }

  a = max((A$mean +A$SD )) *ylimit



  # p=p+Mytheme
  if (cord == TRUE) {
    p<-ggplot(A, aes(x=name , y=mean,group = group ))+
      geom_bar(aes(colour= group,fill = group),stat = "identity", width = 0.9,position = "dodge") +

      geom_errorbar(aes(ymin=mean - SD,
                        ymax=mean + SD,group = group),
                    colour="black",width=0.1,size = 1,position = position_dodge(.9))+

      scale_y_continuous(expand = c(0,0),limits = c(0,a))+#,limits = c(0,a)
      labs(
        y="count"
        # y=name_i
        # title = paste("Normality test",p1,"Homogeneity of variance",p2,sep = ":")
      ) +
      geom_text(data=A, aes(x=name, y=mean +SD ,label=groups),
                position = position_dodge(.9),vjust = -1)+
      guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))
    p
  } else if(cord == FALSE){
    p<-ggplot(A, aes(x=group , y=mean,group = name ))+
      geom_bar(aes(colour= name,fill = name),stat = "identity", width = 0.9,position = "dodge") +

      geom_errorbar(aes(ymin=mean - SD,
                        ymax=mean + SD,group = name),
                    colour="black",width=0.1,size = 1,position = position_dodge(.9))+

      scale_y_continuous(expand = c(0,0),limits = c(0,a))+#,limits = c(0,a)
      labs(
        y="count"
        # y=name_i
        # title = paste("Normality test",p1,"Homogeneity of variance",p2,sep = ":")
      ) +
      geom_text(data=A, aes(x=group, y=mean +SD ,label=groups),position = position_dodge(.9),vjust = -1)+
      guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))
    p
  }

  if (length(unique(data$group))>3){	p=p+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))}
  p


  return(list(p,table = A))
}



