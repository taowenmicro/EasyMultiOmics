#' Variance test or non-parametric test results visualization, using boxplot
#'
#' @param data a data.frame contain the input data
#' @param  i col index wtich need to test
#' @param sig_show Distinctive display, "abc" or "line"
#' @param result output from aovMcomper or KwWlx. You can also import result calculated from other software (a data frame)
#' @param  ns Logical value, whether to display insignificant marks
#' @examples
#' # data(data_wt)
#' result = KwWlx(data = data_wt, i= 4)
#' PlotresultBox = aovMuiBoxP(data = data_wt, i= 4,sig_show ="abc",result = result[[1]])
#' # utput result
#' p = PlotresultBox[[1]]
#' p
#' @return data frame
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


###----使用方差检验结果和多重比较结果做展示：  箱线图展示
aovMuiBoxP = function(data = data_wt, i= 3,sig_show ="line",result = result,ns = FALSE){
  aa = result
  name_i = colnames(data[i])
  data_box <- data %>%
    dplyr::select(1, 2, all_of(i))
  colnames(data_box) = c("ID" , "group","dd" )



  data_box$stat=aa[as.character(data_box$group),]$groups


  max=max(data_box[,c("dd")],na.rm = TRUE)
  min=min(data_box[,c("dd")],na.rm = TRUE)
  x = data_box[,c("group","dd")]
  # y = x %>% group_by(group) %>% summarise(Max=paste('max(',"dd",",na.rm = TRUE",')',sep=""))
  # y=as.data.frame(y)
  # y
  # rownames(y)=y$group
  y = x %>%
    group_by(group) %>%
    summarise(Max = max(dd, na.rm = TRUE))

  y = as.data.frame(y)
  rownames(y) = y$group

  data_box$y=y[as.character(data_box$group),]$Max + (max-min)*0.05

  head(data_box)
  p = ggplot(data_box, aes(x=group, y=.data[["dd"]], color=group)) +
    geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
    labs(
      y=name_i)+
    geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+theme(legend.position="none")



  if (sig_show == "abc") {
    tab.abc = data_box %>% distinct( group, .keep_all = TRUE)
    p = p +
      geom_text(data=tab.abc, aes(x=group, y=y, color=group, label= stat))

    p
  }
  wtq = levels(data$group)
  lis = combn(levels(as.factor(data$group)), 2)
  x <-lis
  my_comparisons <- tapply(x,rep(1:ncol(x),each=nrow(x)),function(i)i)

  line = list()
  if (sig_show == "line") {
    zuhe = combn(aa$group,2)
    xxxx <- tapply(zuhe,rep(1:ncol(zuhe),each=nrow(zuhe)),function(i)i)
    xxxx
    sig_lis = rep("a",dim(zuhe)[2])
    for (i in 1:dim(zuhe)[2]) {
      if (dplyr::filter(aa, group == xxxx[[i]][1])$groups == dplyr::filter(aa, group == xxxx[[i]][2])$groups) {
        sig_lis[i] = "no_sig"
      }

      if (dplyr::filter(aa, group == xxxx[[i]][1])$groups != dplyr::filter(aa, group == xxxx[[i]][2])$groups) {
        sig_lis[i] = "*"
      }

    }

    if (ns == TRUE) {
      #-remove the ns
      xxxx[as.character((1:length(sig_lis))[sig_lis =="no_sig"])] = NULL
      sig_lis = sig_lis[sig_lis != "no_sig"]
    }

    line = list(comparisons = xxxx,annotations=sig_lis,y_position = (seq(from=1, to=max(data_box$dd)/4,length.out=dim(zuhe)[2]) + max(data_box$dd)),tip_length = rep(0.03,dim(zuhe)[2]))
    p = p +
      ggsignif::geom_signif(comparisons = xxxx, annotations=sig_lis,
                  y_position = (seq(from=1, to=max(data_box$dd)/4,length.out=dim(zuhe)[2]) + max(data_box$dd)), tip_length = rep(0.03,dim(zuhe)[2]),color = "black")
    p
  }

  # p=p+Mytheme
  p

  return(list(p,data_box))
}
