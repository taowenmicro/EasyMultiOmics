#' @title Generate a Taxonomic Tree Map for Microbial Composition
#'
#' @description
#' This function constructs a taxonomic tree map from a `phyloseq` object, visualizing microbial taxonomy at various levels (Kingdom to OTU).
#' It uses hierarchical relationships between taxonomic levels to generate a tree structure and calculate abundance values for visualization.
#'
#' @param ps A `phyloseq` object containing microbiome data.
#' @param Top An integer specifying the number of top OTUs to include based on abundance. Default is `200`.
#' @param labtab Optional. A table for labeling the nodes in the taxonomic tree. If `NULL`, default labels are used.
#' @param seed An integer for setting the random seed for reproducibility. Default is `1`.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{A plot object representing the taxonomic tree map with abundance values.}
#'   \item{A data frame containing taxonomic information and abundance values used for plotting.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Filters the input `phyloseq` object to retain the top `Top` OTUs based on abundance.
#'   \item Constructs hierarchical taxonomic labels for each OTU, including all levels from Kingdom to Species.
#'   \item Builds a tree structure using parent-child relationships between taxonomic levels.
#'   \item Calculates mean abundance values for each OTU.
#'   \item Merges taxonomic and abundance data to prepare a data frame for plotting.
#'   \item Generates a taxonomic tree map using either default or custom labels (`labtab`).
#' }
#'
#' This function is particularly useful for exploring microbial community structures and visualizing hierarchical taxonomic relationships.
#'
#' @examples
#' \dontrun{
#' tax = ps.16s %>% vegan_tax() %>%as.data.frame()
#' ps.bac <- ps.16s %>% subset_taxa.wt("Kingdom", "Bacteria")
#' res =maptree.micro (ps = ps.bac,Top = 100,labtab =  NULL,seed = 11)
#' p20 = res[[1]]
#' dat = res[[2]]
#' }
#'
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#'
#' @export

maptree.micro=function (ps = ps, Top = 200, labtab = NULL, seed = 1)
{
  set.seed(seed)
  ps_sub <- filter_OTU_ps(ps = ps, Top = Top)
  ps_sub
  tax = as.data.frame(vegan_tax(ps_sub))
  str(tax)
  tax$Kingdom = paste("D_", tax$Kingdom, sep = "--")
  tax$Phylum = paste(tax$Kingdom, "P_", tax$Phylum, sep = "--")
  tax$Class = paste(tax$Phylum, "C_", tax$Class, sep = "--")
  tax$Order = paste(tax$Class, "O_", tax$Order, sep = "--")
  tax$Family = paste(tax$Order, "F_", tax$Family, sep = "--")
  tax$Genus = paste(tax$Family, "G_", tax$Genus, sep = "--")
  tax$Species = paste(tax$Genus, "S_", tax$Species, sep = "--")
  tax$OTU = paste("ASV", row.names(tax), sep = "--")
  tax_table(ps_sub) = as.matrix(tax)
  head(tax)
  row.names(tax) = tax$OTU
  Desep_group <- as.character(colnames(tax))
  edge_t = tax[c(Desep_group[1:2])]
  dim(edge_t)
  edge_t <- dplyr::distinct(edge_t)
  head(edge_t)
  for (i in 2:7) {
    result = tax[c(Desep_group[i:(i + 1)])]
    colnames(result) = colnames(edge_t)
    result <- dplyr::distinct(result)
    edge_t = rbind(edge_t, result)
  }
  tail(edge_t)
  colnames(edge_t) = c("from", "to")
  if (length(unique(tax$Kingdom)) > 1) {
    edge_t$from = as.character(edge_t$from)
    edge_t$to = as.character(edge_t$to)
    buc = data.frame(from = c("king_up", "king_up"), to = c("D_Archaea",
                                                            "D_Bacteria"))
    row.names(buc) = c("sp_K", "sp_k1")
    deg = rbind(edge_t, buc)
  }
  deg = edge_t
  tree <- FromDataFrameNetwork(deg)
  vertices_t <- data.frame(name = unique(c(as.character(deg$from),
                                           as.character(deg$to))))
  mylevels <- data.frame(name = tree$Get("name"), level = tree$Get("level"))
  vertices_t <- vertices_t %>% dplyr::left_join(., mylevels,
                                                by = c(name = "name"))
  AA = rep("A", length(vertices_t$name))
  for (i in 1:length(vertices_t$name)) {
    AA[i] = strsplit(basename(as.character(vertices_t$name)),
                     "--")[[i]][length(strsplit(basename(as.character(vertices_t$name)),
                                                "--")[[i]])]
  }
  vertices_t$shortName <- AA
  vertices_t <- vertices_t %>% dplyr::mutate(new_label = ifelse(level ==
                                                                  2, shortName, NA))
  row.names(vertices_t) = vertices_t$name
  tax_table = as.data.frame(vegan_tax(ps_sub))
  otu_table = as.data.frame(t(vegan_otu(ps_sub)))
  head(tax_table)
  otu_table$mean = rowMeans(otu_table)
  inde = merge(otu_table, tax_table, by = "row.names", all = TRUE)
  head(inde)
  data_plot <- inde %>% dplyr::group_by(OTU) %>% dplyr::summarise_if(is.numeric,
                                                                     sum) %>% dplyr::right_join(vertices_t, by = c(OTU = "name")) %>%
    as.data.frame()
  data_plot <- data_plot %>% dplyr::left_join(tax_table, by = "OTU") %>%
    dplyr::distinct(OTU, .keep_all = TRUE)
  asa = data_plot$mean
  data_plot$mean[is.na(asa)] = 0
  row.names(data_plot) = data_plot$OTU
  if (is.null(labtab)) {
    p <- maptree_abun(labtab = labtab, deg = deg, data_plot = data_plot)
  }
  else if (!is.null(labtab)) {
    p <- maptree_lab(labtab = labtab, deg = deg, data_plot = data_plot)
  }
  return(list(p, data_plot))
}
