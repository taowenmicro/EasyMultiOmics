#' Variance test or non-parametric test results visualization, using boxbar
#'
#' @param data a data.frame contain the input data
#' @param  i col index wtich need to test
#' @param sig_show Distinctive display, "abc" or "line"
#' @param result output from aovMcomper or KwWlx. You can also import result calculated from other software (a data frame)
#' @param  ns Logical value, whether to display insignificant marks
#' @examples
#' # data(data_wt)
#' result = KwWlx(data = data_wt, i= 4)
#' PlotresultBox = aovMuiBoxBarP(data = data_wt, i= 3,sig_show ="abc",result = result[[1]])
#' # utput result
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


aovMuiBoxBarP = function(data = data_wt, i= 3,sig_show ="line",result = result,ns = FALSE){

  aa = result
  print(i)
  name_i = colnames(data[i])
  name_i
  data_box = data[c(1,2,i)]
  colnames(data_box) = c("ID" , "group","dd" )

  data_box$stat=aa[as.character(data_box$group),]$groups


  max=max(data_box[,c("dd")],na.rm = TRUE)
  min=min(data_box[,c("dd")],na.rm = TRUE)
  x = data_box[,c("group","dd")]

  y = x %>% group_by(group) %>% summarise_(Max=paste('max(',"dd",",na.rm = TRUE",')',sep=""))



  y=as.data.frame(y)
  y
  rownames(y)=y$group
  data_box$y=y[as.character(data_box$group),]$Max + (max-min)*0.05


iris_groups<- group_by(data_box, group)
databar<- dplyr::summarise(iris_groups, mean(dd,na.rm = TRUE), sd(dd,na.rm = TRUE))
colnames(databar) = c("group","mean","sd")

head(databar)
p = ggplot() +
  geom_bar(data = databar,aes(x= group,y = mean,fill = group),stat = "identity", width = 0.4,position = "dodge",colour= "black") +
  geom_errorbar(data = databar,aes(x= group,y = mean,ymin=mean -sd, ymax=mean +sd),colour="black",width=0.1,size = 1)+
  geom_jitter(aes(x=group, y=dd, fill=group),data_box,position=position_jitter(0.17), size=4, alpha=0.7,pch = 21,color = "black")+
  labs(
    y=name_i)
 # p
# geom_hline(aes(yintercept=mean(data_box$dd)), colour="black", linetype=2) +
# geom_vline(aes(xintercept=0), colour="black", linetype="dashed")
head(data_box)
if (sig_show == "abc") {
  p = p +
    geom_text(data=data_box, aes(x=group, y=y*1.1, label= stat))

  p
}


wtq = levels(data$group)
lis = combn(levels(data$group), 2)
x <-lis
my_comparisons <- tapply(x,rep(1:ncol(x),each=nrow(x)),function(i)i)


if (sig_show == "line") {
  zuhe = combn(aa$group,2)
  xxxx <- tapply(zuhe,rep(1:ncol(zuhe),each=nrow(zuhe)),function(i)i)
  xxxx
  sig_lis = rep("a",dim(zuhe)[2])
  for (i in 1:dim(zuhe)[2]) {
    library(tidyverse)

    if (filter(aa, group == xxxx[[i]][1])$groups == filter(aa, group == xxxx[[i]][2])$groups) {
      sig_lis[i] = "no_sig"
    }

    if (filter(aa, group == xxxx[[i]][1])$groups != filter(aa, group == xxxx[[i]][2])$groups) {
      sig_lis[i] = "*"
    }

  }

  # ns = TRUE
  # sig_lis[3] = "no_sig"
  # sig_lis[4] = "no_sig"

  if (ns == TRUE) {
    #-remove the ns
    xxxx[as.character((1:length(sig_lis))[sig_lis =="no_sig"])] = NULL
    sig_lis = sig_lis[sig_lis != "no_sig"]
  }


  p = p +
    ggsignif::geom_signif(aes(x=group, y=dd, fill=group),data_box,comparisons = xxxx, annotations=sig_lis,
                y_position = (seq(from=1, to=max(data_box$dd)/4,length.out=length(sig_lis)) + max(data_box$dd)), tip_length = rep(0.03,length(sig_lis)),color = "black")
  p
}

# p=p+Mytheme
p
return(list(p,data_box,databar))
}

