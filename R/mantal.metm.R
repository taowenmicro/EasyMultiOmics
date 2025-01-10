#' @title Mantel Test and Visualization for Microbial Communities
#'
#' @description
#' The `mantal.metm` function performs a Mantel test to evaluate the correlation between two distance matrices
#' (e.g., microbial community dissimilarities and environmental distances) based on a specified correlation method.
#' It also groups the data by a specified factor and provides a visualization of the results.
#' The function allows customization of the layout of the plots, such as the number of columns and rows.
#'
#' @param ps A `phyloseq` object containing the microbial community data and metadata.
#' @param method A character string specifying the correlation method to be used in the Mantel test.
#'               Options include `"pearson"`, `"spearman"` (default), and `"kendall"`.
#' @param group A character string indicating the name of the grouping variable in the metadata.
#'              This variable will be used to split and analyze the data.
#' @param ncol An integer specifying the number of columns in the plot layout. Default is 5.
#' @param nrow An integer specifying the number of rows in the plot layout. Default is 2.
#'
#' @return This function returns a list containing the following:
#' - Mantel test results for each group.
#' - A visualization of the Mantel test results grouped by the specified factor.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' # Example usage:
#'result <- mantal.metm(ps = ps.16s,method =  "spearman",group = "Group",ncol = gnum,nrow = 1)
#'data <- result[[1]]
#'data
#'p3_7 <- result[[2]]
#' @export

mantal.metm=function (ps = ps, method = "spearman", group = "Group", ncol = 5,
          nrow = 2)
{
  dist <- ggClusterNet::scale_micro(ps = ps, method = "rela") %>%
    ggClusterNet::vegan_otu() %>% vegan::vegdist(method = "bray") %>%
    as.matrix()
  map = phyloseq::sample_data(ps)
  gru = map[, group][, 1] %>% unlist() %>% as.vector()
  id = combn(unique(gru), 2)
  R_mantel = c()
  p_mantel = c()
  name = c()
  R_pro <- c()
  p_pro <- c()
  plots = list()
  for (i in 1:dim(id)[2]) {
    id_dist <- row.names(map)[gru == id[1, i]]
    dist1 = dist[id_dist, id_dist]
    id_dist <- row.names(map)[gru == id[2, i]]
    id_dist = id_dist[1:nrow(dist1)]
    dist2 = dist[id_dist, id_dist]
    mt <- vegan::mantel(dist1, dist2, method = method)
    R_mantel[i] = mt$statistic
    p_mantel[i] = mt$signif
    name[i] = paste(id[1, i], "_VS_", id[2, i], sep = "")
    mds.s <- vegan::monoMDS(dist1)
    mds.r <- vegan::monoMDS(dist2)
    pro.s.r <- vegan::protest(mds.s, mds.r)
    R_pro[i] <- pro.s.r$ss
    p_pro[i] <- pro.s.r$signif
    Y <- cbind(data.frame(pro.s.r$Yrot), data.frame(pro.s.r$X))
    X <- data.frame(pro.s.r$rotation)
    Y$ID <- rownames(Y)
    p1 <- ggplot(Y) + geom_segment(aes(x = X1, y = X2, xend = (X1 +
                                                                 MDS1)/2, yend = (X2 + MDS2)/2), arrow = arrow(length = unit(0,
                                                                                                                             "cm")), color = "#B2182B", size = 1) + geom_segment(aes(x = (X1 +
                                                                                                                                                                                            MDS1)/2, y = (X2 + MDS2)/2, xend = MDS1, yend = MDS2),
                                                                                                                                                                                 arrow = arrow(length = unit(0, "cm")), color = "#56B4E9",
                                                                                                                                                                                 size = 1) + geom_point(aes(X1, X2), fill = "#B2182B",
                                                                                                                                                                                                        size = 4, shape = 21) + geom_point(aes(MDS1, MDS2),
                                                                                                                                                                                                                                           fill = "#56B4E9", size = 4, shape = 21) + labs(title = paste(id[1,
                                                                                                                                                                                                                                                                                                           i], "-", id[2, i], " ", "Procrustes analysis:\n    M2 = ",
                                                                                                                                                                                                                                                                                                        round(pro.s.r$ss, 3), ", p-value = ", round(pro.s.r$signif,
                                                                                                                                                                                                                                                                                                                                                    3), "\nMantel test:\n    r = ", round(R_mantel[i],
                                                                                                                                                                                                                                                                                                                                                                                          3), ", p-value =, ", round(p_mantel[i], 3), sep = ""))
    p1
    plots[[i]] = p1
  }
  dat = data.frame(name, R_mantel, p_mantel, R_pro, p_pro)
  pp = ggpubr::ggarrange(plotlist = plots, common.legend = TRUE,
                         legend = "right", ncol = ncol, nrow = nrow)
  return(list(dat, pp))
}
