#' Stacked histograms add error bars, marked by significance
#'
#' @param data a data.frame contain the input data
#' @param  i col index wtich need to test
#' @param errbar add a distinctive label, TRUE or FEASE could be selected
#' @param result output from aovMcomper or KwWlx. You can also import result calculated from other software (a data frame)
#' @param add_abc add sig abc T or F
#' @examples
#' # data(data_wt)
#'result = MuiKwWlx(data = data_wt,num = c(4,9,8))
#'res <- value_stackBar(data = data_wt,i = c(4,9,8) ,result = result,
#'                      add_abc = TRUE)
#'res[[1]]
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} PengHao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export



value_stackBar = function(
  data = data,
  i,
  result = result,
  add_abc = TRUE
){
  i = c(1,2,i)
  data <- data[i]
  tem = data[,c(-1)]
  colnames(tem)[1] = "id"


  result$group = row.names(result)
  abc <- reshape2::melt(result, id="group",
                        variable.name="variable", value.name = "abc")
  colnames(abc)[1] = "id"
  head(tem)

  mean = function(x) {
    base::mean(x,na.rm = TRUE)
  }
  tem.2 = tem %>% dplyr::group_by(id) %>%
    dplyr::summarise_if(is.numeric, mean) %>%
    as.data.frame()
  count = tem.2[,c(-1)]
  tem.2[,c(-1)] = t(t(count)/colSums(count))

  pcm = reshape2::melt(tem.2, id = c("id"))

  dat = pcm %>% dplyr::left_join(abc)
  if (add_abc == FALSE) {
    p1 = ggplot(dat, aes(x = value, y = variable,group = id,fill = id)) +
      geom_bar(stat="identity",color="black", width=.6) +
      theme(
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.text.x = element_text(colour = "black",angle = 90)
      )
    tem.6 = dat
  }

  if (add_abc == TRUE) {
    ## Data statistics mean, standard deviation, standard error
    labfun = function(x) {
      cumsum(x) - 0.5*x
    }
    tem.3 <- stats::aggregate(pcm$value, by=list(pcm$variable),
                              FUN= labfun)

    tem.3 = data.frame(group = tem.3$Group.1,tem.3$x)
    tem.id = pcm$id %>% unique() %>% as.character()
    colnames(tem.3) = c("variable",tem.id)
    head(tem.3)
    tem.4 = reshape2::melt(tem.3, id = c("variable"))
    colnames(tem.4) = c("variable","id","label_x")

    tem.5 = pcm %>%
      dplyr::left_join(abc)

    tem.6 = tem.5 %>%
      dplyr::left_join(tem.4)

    head(tem.6)
    tem.6$id = factor(tem.6$id,levels = tem.id[length(tem.id):1] )


    p1 = ggplot(tem.6, aes(x = value, y = variable,group = id,fill = id)) +
      geom_bar(stat="identity",color="black", width=.6) +
      geom_text(aes(x = label_x, label = abc),size = 2) +
      theme(
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.text.x = element_text(colour = "black",angle = 90)
      )
  }

  return(list(p1,tem.6))
}

