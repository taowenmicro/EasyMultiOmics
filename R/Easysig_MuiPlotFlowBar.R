#' Stacked histograms add error bars, marked by significance
#'
#' @param data a data.frame contain the input data
#' @param  i col index wtich need to test
#' @param errbar add a distinctive label, TRUE or FEASE could be selected
#' @param result output from aovMcomper or KwWlx. You can also import result calculated from other software (a data frame)
#' @param add_abc add sig abc T or F
#' @param order order with colnames according to value
#' @param sample T or F show sample data not group data
#' @examples
#' # data(data_wt)
#' result = MuiKwWlx(data = data_wt,num = c(4,9,8))
#' res <- MuiPlotFlowBar(data = data_wt,i = c(4,9,8) ,result = result)
#' # utput result
#' res[[1]]
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} PengHao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export




MuiPlotFlowBar <- function(data = data,
                            i,
                            result = result,
                            errbar = TRUE,
                            add_abc = TRUE,
                            order = FALSE,
                            sample = FALSE,
                           flow = F,
                            ...){

  #-- data prepare
  i = c(2,i)
  data <- data[i]
  if (order == TRUE) {
    tem = colSums(data[,-1]) %>%
      as.data.frame()
    colnames(tem) = "sum"
    id = tem %>%
      dplyr::arrange(sum) %>%
      row.names()
    data= cbind(group = data[,1],data[,-1][,c(match(id,colnames(data[,-1])))])
  }
  if (sample == TRUE) {
    data$group = row.names(data)
    df <- reshape2::melt(data, id="group", variable.name="variable", value.name = "Size")
    colnames(df) = c("Sample","Rank","Relative_abundance")
    head(df)
    p <-  ggplot(df, aes(x= Sample, y= Relative_abundance, fill=Rank)) +
      geom_bar(stat="identity",color="black", width=.6)
    plotdata = df



  } else {
    result$group = row.names(result)
    abc <- reshape2::melt(result, id="group", variable.name="variable", value.name = "abc")

    # data for plot
    df <- reshape2::melt(data, id="group", variable.name="variable", value.name = "Size")
    head(df)
    ## Data statistics mean, standard deviation, standard error
    mean <- stats::aggregate(df$Size, by=list(df$group, df$variable), FUN=mean)
    sd <- stats::aggregate(df$Size, by=list(df$group, df$variable), FUN=sd)
    len <- stats::aggregate(df$Size, by=list(df$group, df$variable), FUN=length)
    df_res <- data.frame(mean, sd=sd$x, len=len$x)
    colnames(df_res) = c("group", "variable", "Mean", "Sd", "Count")
    df_res$Se <- df_res$Sd/sqrt(df_res$Count)
    levels(df_res $variable) = as.character(unique(df_res$variable))

    # Construct error line coordinates--
    # df_res = plyr::ddply(df_res,"group",transform,label_y = cumsum(Mean))
    # Construct distinctive marker positions
    df_res_sub = plyr::ddply(df_res,"group", summarize,label_y = cumsum(Mean),  label_abc = cumsum(Mean) - 0.5*Mean,
                             variable = variable)
    # df_res = cbind(df_res,df_res_sub[-1])
    df_res <- df_res%>%
      dplyr::inner_join(df_res_sub )
    # Factor rearrangement
    df_res$variable = factor(df_res$variable,order = F,
                             levels = levels(df_res$variable)[length(levels(df_res$variable)):1])
    #--conbind plot data
    plotdata <- df_res %>%
      dplyr::left_join(abc,by = c("group","variable"))

    plotdata$variable = factor(plotdata$variable,
                               levels = as.character(unique(df_res$variable))[length(levels(df_res$variable)):1])

    cs = plotdata$variable
    lengthfactor <- plotdata$variable %>%
      levels() %>%
      length()
    cs4 <- cs %>%
      as.factor() %>%
      summary()  %>%
      as.data.frame()
    cs4$id = row.names(cs4)


    df_arrange<- dplyr::arrange(cs4, id)
    Taxonomies_x1<- dplyr::arrange(plotdata , variable)
    Taxonomies_x1$ID = factor(rep(c(1:lengthfactor), cs4$.))

    if (flow) {
      # flower plot
      p <- ggplot(Taxonomies_x1,
                   aes(x = group, alluvium = variable, y = Mean)) +
        ggalluvial::geom_flow(aes(fill = variable, colour = variable), width = 0)

    }else {
      p <- ggplot(Taxonomies_x1, aes(x = group, y = Mean,fill = variable,
                                     alluvium = variable,stratum = ID)) +
        ggalluvial::geom_flow(aes(fill = variable, colour = aa),
                              stat = "alluvium", lode.guidance = "rightleft",
                              color = "black",size = 0.2,width = 0.35,alpha = .2)  +
        geom_bar(width = 0.45,stat = "identity")
    }

    if (add_abc == TRUE) {
      p = p +
        geom_text(aes(y = label_abc, label = abc))
    }


    if (errbar == TRUE) {
      p <- p + geom_errorbar(aes(ymin=label_y-Sd, ymax=label_y +Sd), width=.2)
    }

  }
  # id = plotdata %>%
  #   dplyr::group_by(variable) %>%
  #   dplyr::summarise(mean(Mean)) %>%
  #   dplyr::arrange(`mean(Mean)`) %>%
  #   .$variable
  # plotdata$variable = factor(as.character(plotdata$variable),
  #                            levels = as.character(id))


  return(list(p,plotdata))
}

