#' Visualize multiple sets of data Scatter plot + barplot
#'
#' @param data a data.frame contain the input data
#' @param  num col index wtich need to test
#' @param sig_show Distinctive display, "abc" or "line"
#' @param result output from aovMcomper or KwWlx. You can also import result calculated from other software (a data frame)
#' @param path creat a folder to save the single plot.
#' @examples
#' # data(data_wt)
#' result = KwWlx(data = data_wt, i= 4)
#' MuiPlotReBoxBar(data = data_wt,num = c(4:6),result = result,sig_show ="abc")
#'
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export

MuiPlotReBoxBar = function(data = data_wt,num = c(4:6),result = result,sig_show ="abc",path = "./MuiboxBar/"){
  dir.create(path)
  data = data
  Mytheme <- theme_bw()+

    # scale_fill_manual(values = mi, guide = guide_legend(title = NULL))+
    theme(

      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),

      plot.title = element_text(vjust = -8.5,hjust = 0.1),
      axis.title.y =element_text(size = 20,face = "bold",colour = "black"),
      axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
      axis.text = element_text(size = 20,face = "bold"),
      axis.text.x = element_text(colour = "black",size = 14),
      axis.text.y = element_text(colour = "black",size = 14),
      legend.text = element_text(size = 15,face = "bold"),
      legend.position = "none"#是否删除图例

    )
  # N = 3

  for (N in num) {
    name = colnames(data[N])

    as = result[match( name,colnames(result))]
    as
    colnames(as) = "groups"
    as$group = row.names(as)

    PlotresultBox = aovMuiBoxBarP(data = data, i= N,sig_show =sig_show,result = as)
    p = PlotresultBox[[1]]
    p
    if (dim(as)[1]>3){	p=p+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))}
    name = colnames(data[N])

    path = path
    dir.create(path)
    FileName <- paste(path,name,"_boxbar", ".pdf", sep = "_")
    ggsave(FileName, p, width = 8, height = 8)

  }
}
