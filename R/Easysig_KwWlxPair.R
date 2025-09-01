#'  Easy for Non-parametric test
#'
#' @param data a data.frame contain the input data
#' @param  i col index wtich need to test
#' @param  method "t.test" or "wilcox.test" could be selected by method
#' @param paired a logical indicating whether you want a paired test. Used in t.test and in wilcox.test.
#' @examples
#' # data(data_wt)
#' result = KwWlx2(data = data_wt, i= 4,method = "t.test",paired = FALSE)
#' result[[1]]
#' @return data frame
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{2018203048@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


#-------------------------------------------Non-parametric test-------------------------------
KwWlxPair = function(data = data_wt, i= 3,method = "wilcox.test",paired = FALSE){
  # i = 17
  ss <- data %>%
    dplyr::select("group",count = i)

  # kruskal.test
  krusk = ggpubr::compare_means(count ~ group, data=ss, method = "kruskal.test",paired = paired)
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
  aa
  # aa =ord_sig(data = aa,ID = "groups")
  dat <- ss %>% group_by(group) %>%
    summarise(mean = mean(count)) %>%
    inner_join(as_tibble(aa),by = c("group" = "group")) %>%
    arrange(desc(mean))
  tmp.1 <- dat$groups %>% unique()
  tmp.2 <- data.frame(ori = sort(tmp.1,decreasing =TRUE),new = sort(tmp.1,decreasing =FALSE))
  i = 1
  for (i in 1:nrow(dat)) {
    dat[i,3] <- tmp.2$new[match(dat[i,3],tmp.2$ori)]
  }
  dat$mean = NULL
  dat = dat[match(dat$group,aa$group),][,c(2,1)] %>% as.data.frame()
  row.names(dat) = dat$group

  return(list(aa,wilcox = krusk,kruskal = sumkrusk))
}
