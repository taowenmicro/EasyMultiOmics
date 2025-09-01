#' Variance test or non-parametric test results visualization, using histogram
#'
#' @param data a data.frame contain the input data
#' @param  i col index wtich need to test
#' @param sig_show Distinctive display, "abc" or "line"
#' @param result output from aovMcomper or KwWlx. You can also import result calculated from other software (a data frame)
#' @param  ns Logical value, whether to display insignificant marks
#' @examples
#' # data(data_wt)
#' result = KwWlx(data = data_wt, i= 4)
#' PlotresultBar = aovMuiBarPlot(data = data_wt, i= 4,sig_show ="abc",result = result[[1]])
#' # utput result
#' PlotresultBar[[1]]
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export



###----使用方差检验结果和多重比较结果做展示：  柱状图展示
aovMuiBarPlot = function(data = data_wt, i= 3,sig_show ="line",result = result,ns = FALSE){

  name_i = colnames(data[i])
  # the mean and variance
  wen1 = as.data.frame(tapply(as.vector(as.matrix(data[i])),data$group,mean,na.rm=TRUE))
  wen2 = as.data.frame(tapply(as.vector(as.matrix(data[i])),data$group,sd,na.rm=TRUE))
  went = cbind(wen1,wen2)
  colnames(went) = c("mean" ,"SD")
  aa = result
  wentao = merge(aa,went, by="row.names",all=F)

  aa = mutate(wentao, ymin = mean - SD, ymax =  mean + SD)
  a = max(aa$ymax)*1.2##

  ###
  p = ggplot(aa , aes(x = group, y = mean,colour= group)) +
    geom_bar(aes(colour= group,fill = group),stat = "identity", width = 0.4,position = "dodge") +

    geom_errorbar(aes(ymin=ymin,
                      ymax=ymax),
                  colour="black",width=0.1,size = 1)+
    scale_y_continuous(expand = c(0,0),limits = c(0,a))+#
    labs(
      # x=paste(name_i,"of all group", sep = "_"),
      y=name_i
      # title = paste("Normality test",p1,"Homogeneity of variance",p2,sep = ":")
    )
  line = list()
  if (sig_show == "line") {
    zuhe = combn(aa$group,2)
    xxxx <- tapply(zuhe,rep(1:ncol(zuhe),each=nrow(zuhe)),function(i)i)
    xxxx
    sig_lis = rep("a",dim(zuhe)[2])
    for (i in 1:dim(zuhe)[2]) {


      if (filter(aa, group == xxxx[[i]][1])$groups == filter(aa, group == xxxx[[i]][2])$groups) {
        sig_lis[i] = "no_sig"
      }

      if (filter(aa, group == xxxx[[i]][1])$groups != filter(aa, group == xxxx[[i]][2])$groups) {
        sig_lis[i] = "*"
      }

    }

    if (ns == TRUE) {
      #-remove the ns
      xxxx[as.character((1:length(sig_lis))[sig_lis =="no_sig"])] = NULL
      sig_lis = sig_lis[sig_lis != "no_sig"]
    }

    line = list(comparisons = xxxx,annotations=sig_lis,y_position = (seq(from=1, to=max(aa$mean)/4,length.out=dim(zuhe)[2]) + max(aa$mean)),tip_length = rep(0.03,dim(zuhe)[2]))

    p = p +
      ggsignif::geom_signif(comparisons = xxxx, annotations=sig_lis,
                  y_position = (seq(from=1, to=max(aa$mean)/4,length.out=dim(zuhe)[2]) + max(aa$mean)), tip_length = rep(0.03,dim(zuhe)[2]),color = "black")
    p
  }


  if (sig_show == "abc") {


    p = p + geom_text(aes(label = groups,y=ymax, x = group,vjust = -0.3,size = 6))
    p
  }

  p

  if (length(unique(data$group))>3){	p=p+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))}

  return(list(p,wentao,aa,line))
}
