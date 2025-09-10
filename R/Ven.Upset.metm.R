#' @title Generate Venn Diagrams and Upset Plots for Microbial Group Comparisons
#' @description
#' This function creates Venn diagrams and Upset plots to compare microbial groups based on presence/absence thresholds.
#'
#' @param otu An OTU table (optional). If `NULL`, uses the `ps` object.
#' @param tax A taxonomy table (optional). If `NULL`, uses the `ps` object.
#' @param map A sample metadata table (optional). If `NULL`, uses the `ps` object.
#' @param tree A phylogenetic tree (optional). If `NULL`, uses the `ps` object.
#' @param ps A `phyloseq` object containing microbiome data.
#' @param group A character string specifying the grouping variable in the sample metadata. Default is `"Group"`.
#' @param N A numeric threshold for presence-absence data. Default is `0.5`.
#' @param size Numeric. The scaling factor for the Upset plot. Default is `3`.
#'
#' @return A list containing:
#' \describe{
#'   \item{p}{A `ggVennDiagram` object representing the Venn diagram.}
#'   \item{p2}{A composite Upset plot.}
#'   \item{ven2}{A processed binary presence-absence data frame.}
#'   \item{main_plot}{The main bar plot of the Upset plot.}
#'   \item{side_plot}{The side bar plot of the Upset plot.}
#' }
#' @author
#' Tao  Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @export
#' @examples
#' \dontrun{
#' library(phyloseq)
#' # Example with a phyloseq object
#' res = Ven.Upset.metm(ps =  ps.16s,group = "Group",N = 0.5,size = 3)
#' p10.1 = res[[1]]
#' p10.1
#' p10.2 = res[[2]]
#' p10.2
#' p10.3 = res[[4]]
#' p10.3
#' dat = res[[3]]
#' head(dat)
#' }
#'



Ven.Upset.metm =function (otu = NULL, tax = NULL, map = NULL, tree = NULL, ps = NULL,
                          group = "Group", N = 0.5, size = 3)
{
  library(ggVennDiagram)
  ps = ggClusterNet::inputMicro(otu, tax, map, tree, ps, group = group)
  aa = ggClusterNet::vegan_otu(ps)
  otu_table = as.data.frame(t(aa))
  count = aa
  countA = count
  sub_design <- as.data.frame(phyloseq::sample_data(ps))
  count[count > 0] <- 1
  count2 = as.data.frame(count)
  aa = sub_design[, "Group"]
  colnames(aa) = "Vengroup"
  iris.split <- split(count2, as.factor(aa$Vengroup))
  iris.apply <- lapply(iris.split, function(x) colSums(x[]))
  iris.combine <- do.call(rbind, iris.apply)
  ven2 = t(iris.combine)
  for (i in 1:length(unique(phyloseq::sample_data(ps)[, "Group"]))) {
    aa <- as.data.frame(table(phyloseq::sample_data(ps)[,
                                                        "Group"]))[i, 1]
    bb = as.data.frame(table(phyloseq::sample_data(ps)[,
                                                       "Group"]))[i, 2]
    ven2[, aa] = ven2[, aa]/bb
  }
  ven2[ven2 < N] = 0
  ven2[ven2 >= N] = 1
  ven2 = as.data.frame(ven2)
  ven3 = list()
  for (i in 1:ncol(ven2)) {
    ven3[[i]] = row.names(ven2)[ven2[, i] == 1]
  }
  names(ven3) = colnames(ven2)
  p = ggVennDiagram(ven3, label_alpha = 0, category.names = names(ven3)) +
    ggplot2::scale_fill_gradientn(colours = c("#E41A1C",
                                              "#377EB8", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
                                              "#999999"))
  A = list()
  for (i in 1:nrow(ven2)) {
    tem = ven2[i, ] %>% t() %>% .[, 1]
    A[[i]] = tem[tem == 1] %>% names()
  }
  dat = tibble(ID = row.names(ven2), A)
  head(dat)
  main_plot <- ggplot(data = dat, aes(x = A)) + geom_bar() +
    ggupset::scale_x_upset(n_intersections = 20) + theme_classic()
  side_plot <- dat %>% dplyr::select(A) %>% unnest(cols = A) %>% dplyr::count(A) %>%
    mutate(A = fct_reorder(as.factor(A), n)) %>% ggplot(aes(y = n,
                                                            x = A)) + geom_col() + coord_flip() + scale_y_reverse() +
    xlab("") + ylab("") + theme(axis.ticks.y = element_blank(),
                                axis.text.y = element_blank(), panel.background = element_blank(),
                                panel.grid = element_blank())
  p2 = cowplot::plot_grid(cowplot::plot_grid(NULL, side_plot +
                                               theme(plot.margin = unit(c(1, -5, -5, 1), "pt")), ncol = 1,
                                             rel_heights = c(size, 1)), main_plot, nrow = 1, rel_widths = c(1,
                                                                                                            3))
  return(list(p, p2, ven2, main_plot, side_plot))
}
