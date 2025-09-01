

FacetMuiPlotresultBar2 = function (data = data_wt, num = c(4:6), result = result, sig_show = "abc", 
          ncol = 3, fac.level = NULL) 
{
  N = num[1]
  name = colnames(data[N])
  as = result[match(name, colnames(result))]
  colnames(as) = "groups"
  as$group = row.names(as)
  PlotresultBar = aovMuiBarPlot(data = data, i = N, sig_show = sig_show, 
                                result = as)
  p = PlotresultBar[[2]]
  p
  name = colnames(data[N])
  p$name = name
  A = p
  for (N in num[-1]) {
    name = colnames(data[N])
    as = result[match(name, colnames(result))]
    colnames(as) = "groups"
    as$group = row.names(as)
    PlotresultBox = aovMuiBarPlot(data = data, i = N, sig_show = sig_show, 
                                  result = as)
    p = PlotresultBox[[2]]
    p
    name = colnames(data[N])
    p$name = name
    A = rbind(A, p)
  }
  if (!is.null(fac.level)) {
    A$name = factor(A$name, levels = fac.level)
  }
  p <- ggplot(A, aes(y = group, x = mean)) + geom_bar(aes(fill = group), 
                                                      stat = "identity",
                                                      width = 0.4, position = "dodge") + 
    geom_bar(data = A, aes(y = 1, x = (mean + SD) * 1.1), 
             stat = "identity", width = 0.4, position = "dodge", 
             alpha = 0) + geom_errorbar(aes(xmin = mean - SD, 
                                            xmax = mean + SD), colour = "black",
                                        width = 0.1, 
                                        size = 1) + 
    scale_x_continuous(expand = c(0, 0)) + labs() + 
    geom_text(data = A, aes(y = group, x = mean + SD, label = groups), 
              hjust = -1) + guides(color = guide_legend(title = NULL), 
                                   shape = guide_legend(title = NULL)) +
    facet_wrap(. ~ name, scales = "free_x", ncol = ncol)
  p
  if (length(unique(data$group)) > 3) {
    p = p + theme(axis.text.y = element_text(angle = 45, 
                                             vjust = 1, hjust = 1))
  }
  p
  return(list(p, table = A))
}
