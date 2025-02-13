#' CNPS Gene Network Analysis for Phyloseq Object
#'
#' This function processes a `phyloseq` object by mapping the OTU data to CNPS gene groups,
#' performs correlation-based network analysis, and visualizes the network using ggplot2.
#' It returns a list containing the network plot, nodes, and edges.
#'
#' @param ps A `phyloseq` object (default is `ps.kegg`) containing OTU data to be mapped to CNPS gene groups.
#' @param id.0 A string representing the CNPS gene group to focus on (default is `"C"`).
#'
#' @return A list containing:
#'   \itemize{
#'     \item A ggplot object (`p1`) representing the CNPS gene network plot.
#'     \item A data frame (`nodes`) containing the nodes of the network.
#'     \item A data frame (`edge`) containing the edges of the network.
#'   }
#' @details
#' This function performs the following steps:
#' \itemize{
#'   \item Filters the CNPS gene data to focus on a specific gene group (`id.0`).
#'   \item Extracts OTU data from the `phyloseq` object and transforms it.
#'   \item Subsets the OTU data for the specified CNPS gene group.
#'   \item Performs correlation-based analysis on the OTU data using Pearson's method.
#'   \item Creates a network graph using the `igraph` package based on the correlation matrix.
#'   \item Identifies and visualizes the network using ggplot2, with nodes representing OTUs and edges representing correlations.
#' }
#' @examples
#' # Assuming you have a phyloseq object `ps.kegg`
#' network_results <-CNPS.network2(ps = ps.kegg,id.0 = "C")
#' network_plot <- network_results[[1]]
#'
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
CNPS.network2 = function (ps = ps.kegg, id.0 = "C")
{
  dat = EasyMultiOmics.db::db.cnps
  id = dat %>% filter(Group == id.0) %>% .$K.number
  tax = dat %>% filter(Group == id.0) %>% distinct(K.number,
                                                   .keep_all = TRUE)
  head(tax)
  row.names(tax) = tax$K.number
  tax$symbol %>% table()
  ps.t = ps %>% scale_micro(method = "rela") %>% subset_taxa.wt2("OTU",
                                                                 id)
  tax_table(ps.t) = as.matrix(tax)
  ps.t2 = ps.t %>% tax_glom_wt(ranks = "symbol")
  map = sample_data(ps.t2)
  head(map)
  map.t = ps.t2 %>% sample_data() %>% as.data.frame()
  gru3 = map.t$Group %>% unique() %>% as.character()
  n = 2
  for (n in 1:length(gru3)) {
    ps.t4 = ps.t2 %>% subset_samples.wt("Group", gru3[n]) %>%
      filter_taxa(function(x) sum(x) > 0, TRUE)
    result = cor_Big_micro(ps = ps.t4, N = 0, r.threshold = 0.6,
                           p.threshold = 0.05, method = "pearson")
    cor = result[[1]]
    dim(cor)
    igraph = make_igraph(cor)
    library(igraph)
    dat = net_properties.4(igraph, n.hub = F)
    head(dat, n = 16)
    netClu = data.frame(ID = row.names(cor), group = rep(1,
                                                         length(row.names(cor)))[1:length(row.names(cor))])
    netClu$group = as.factor(netClu$group)
    result2 = PolygonClusterG(cor = cor, nodeGroup = netClu)
    node = result2[[1]]
    edge = edgeBuild(cor = cor, node = node)
    head(edge)
    if (dim(edge)[1] != 0) {
      edge2 = edge %>% filter(weight != 0)
      tem = c(edge2$OTU_2, edge2$OTU_1) %>% unique()
      cor2 = cor[tem, tem]
      netClu = data.frame(ID = row.names(cor2), group = rep(1,
                                                            length(row.names(cor2)))[1:length(row.names(cor2))])
      netClu$group = as.factor(netClu$group)
      result2 = PolygonClusterG(cor = cor2, nodeGroup = netClu)
      node = result2[[1]]
      otu = ps.t4 %>% vegan_otu() %>% t() %>% as.data.frame()
      head(otu)
      otu$mean = rowSums(otu)
      otu = otu %>% rownames_to_column("ID")
      nodes = node %>% left_join(otu, by = c(elements = "ID"))
      head(nodes)
      edge = edgeBuild(cor = cor2, node = node)
      head(edge)
      edge$abs = abs(edge$weight)
      library(ggrepel)
      library(ggnewscale)
      p1 <- ggplot() + geom_segment(aes(x = X1, y = Y1,
                                        xend = X2, yend = Y2, color = abs), data = edge,
                                    size = 1) + geom_point(aes(X1, X2, size = mean),
                                                           pch = 21, data = nodes, fill = "#984EA3") + scale_color_gradientn(colours = c("white",
                                                                                                                                         "grey10")) + scale_x_continuous(breaks = NULL) +
        scale_y_continuous(breaks = NULL) + ggrepel::geom_text_repel(aes(X1,
                                                                         X2, label = elements), size = 4, data = nodes) +
        theme(panel.background = element_blank()) + theme(axis.title.x = element_blank(),
                                                          axis.title.y = element_blank()) + theme(legend.background = element_rect(colour = NA)) +
        theme(panel.background = element_rect(fill = "white",
                                              colour = NA)) + theme(panel.grid.minor = element_blank(),
                                                                    panel.grid.major = element_blank())
      p1
    }
  }
  return(list(p1, nodes, edge))
}
