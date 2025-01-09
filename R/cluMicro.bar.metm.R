
#' @title Generate Clustered Bar Plots for Microbial Data
#'
#' @description
#' This function performs hierarchical clustering on microbial data and generates stacked bar plots and dendrograms to visualize the clustering results. The microbial data can be grouped at a specified taxonomic rank and converted to relative abundances for visualization.
#'
#' @param dist A character string specifying the distance method for clustering. Default is `"bray"`.
#' @param otu A data frame containing OTU data. If `NULL`, the `ps` object is used.
#' @param tax A data frame containing taxonomy data. If `NULL`, the `ps` object is used.
#' @param map A data frame containing sample metadata. If `NULL`, the `ps` object is used.
#' @param tree A phylogenetic tree object. If `NULL`, the `ps` object is used.
#' @param j A character string specifying the taxonomic rank to analyze (e.g., `"Phylum"`, `"Genus"`). Default is `"Phylum"`.
#' @param ps A `phyloseq` object containing microbiome data.
#' @param rep An integer specifying the number of replicates. Default is `6`.
#' @param Top An integer specifying the number of top taxa to display. Taxa not in the top are grouped as `"Other"`. Default is `10`.
#' @param tran Logical. If `TRUE`, transforms data to relative abundances. Default is `TRUE`.
#' @param hcluter_method A character string specifying the hierarchical clustering method. Default is `"complete"`.
#' @param Group A character string specifying the grouping variable in the sample metadata. Default is `"Group"`.
#' @param cuttree An integer specifying the number of clusters to cut the dendrogram into. Default is `3`.
#'
#' @return
#' A list containing the following elements:
#' \describe{
#'   \item{A `ggtree` dendrogram showing the clustering tree at the sample level.}
#'   \item{A stacked bar plot combined with the dendrogram at the sample level.}
#'   \item{A `ggtree` dendrogram showing the clustering tree at the group level.}
#'   \item{A stacked bar plot combined with the dendrogram at the group level.}
#'   \item{A data frame with microbial composition data, including relative abundances and taxonomic annotations.}
#' }
#'
#' @details
#' This function:
#' \itemize{
#'   \item Standardizes the microbial data to relative abundances if `tran = TRUE`.
#'   \item Performs hierarchical clustering on both sample-level and group-level data using the specified distance and clustering methods.
#'   \item Combines the clustering dendrograms with stacked bar plots of microbial composition at the specified taxonomic rank (`j`).
#'   \item Groups microbial taxa not in the top `Top` by abundance as `"Other"`.
#' }
#'
#' The function outputs plots and processed data to explore microbial community composition and clustering relationships.
#'
#' @examples
#' \dontrun{result <-  cluMicro.bar.metm (dist = "bray",ps = ps.16s,j = "Genus",Top = 10, # 提取丰度前十的物种注释tran = TRUE, # 转化为相对丰度值hcluter_method = "complete",Group = "Group",cuttree = length(unique(phyloseq::sample_data(ps)$Group)))
#' result[[1]]
#' p5_2 <- result[[2]]
#' p5_3 <- result[[3]]
#' p5_4 <- result[[4]]
#' clubardata <- result[[5]]
#' }
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @export


cluMicro.bar.metm=function (dist = "bray", otu = NULL, tax = NULL, map = NULL,
          tree = NULL, j = "Phylum", ps = ps, rep = 6, Top = 10, tran = TRUE,
          hcluter_method = "complete", Group = "Group", cuttree = 3)
{
  ps = ggClusterNet::inputMicro(otu, tax, map, tree, ps, group = Group)
  ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x/sum(x))
  otu = as.data.frame(t(ggClusterNet::vegan_otu(ps1_rela)))
  unif = phyloseq::distance(ps1_rela, method = dist)
  hc <- stats::hclust(unif, method = hcluter_method)
  clus <- stats::cutree(hc, cuttree)
  d = data.frame(label = names(clus), member = factor(clus))
  map = as.data.frame(phyloseq::sample_data(ps))
  dd = merge(d, map, by = "row.names", all = F)
  row.names(dd) = dd$Row.names
  dd$Row.names = NULL
  p = ggtree::ggtree(hc) %<+% dd + geom_tippoint(size = 5,
                                                 shape = 21, aes(fill = Group, x = x)) + geom_tiplab(aes(color = Group,
                                                                                                         x = x * 1.2), hjust = 1)
  p
  psdata = ggClusterNet::tax_glom_wt(ps = ps1_rela, ranks = j)
  if (tran == TRUE) {
    psdata = psdata %>% phyloseq::transform_sample_counts(function(x) {
      x/sum(x)
    })
  }
  otu = phyloseq::otu_table(psdata)
  tax = phyloseq::tax_table(psdata)
  for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {
      tax[i, j] = tax[i, j]
    }
    else {
      tax[i, j] = "Other"
    }
  }
  phyloseq::tax_table(psdata) = tax
  Taxonomies <- psdata %>% phyloseq::psmelt()
  head(Taxonomies)
  Taxonomies$Abundance = Taxonomies$Abundance * 100
  Taxonomies$OTU = NULL
  colnames(Taxonomies)[1] = "id"
  head(Taxonomies)
  p <- p + ggnewscale::new_scale_fill()
  p
  p1 <- facet_plot(p, panel = "Stacked Barplot", data = Taxonomies,
                   geom = ggstance::geom_barh, mapping = aes(x = Abundance,
                                                             fill = !!sym(j)), color = "black", stat = "identity")
  p1
  grotax <- Taxonomies %>% dplyr::group_by(Group, !!sym(j)) %>%
    dplyr::summarise(Abundance = sum(Abundance))
  head(grotax)
  data = c()
  i = 2
  for (i in 1:length(unique(phyloseq::sample_data(psdata)$Group))) {
    a <- as.data.frame(table(phyloseq::sample_data(psdata)$Group))[i,
                                                                   1]
    b = as.data.frame(table(phyloseq::sample_data(psdata)$Group))[i,
                                                                  2]
    c <- grotax %>% dplyr::filter(Group == a)
    c$Abundance <- c$Abundance/b
    head(c)
    data = c
    if (i == 1) {
      table = data
    }
    if (i != 1) {
      table = rbind(table, data)
    }
  }
  sum(grotax$Abundance)
  head(table)
  ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x/sum(x))
  otu = as.data.frame((ggClusterNet::vegan_otu(ps1_rela)))
  iris.split <- split(otu, as.factor(phyloseq::sample_data(ps1_rela)$Group))
  iris.apply <- lapply(iris.split, function(x) colMeans(x))
  iris.combine <- do.call(rbind, iris.apply)
  otuG = t(iris.combine)
  ps = phyloseq::phyloseq(phyloseq::otu_table(otuG, taxa_are_rows = T),
                          phyloseq::tax_table(ps1_rela))
  hc = ps %>% phyloseq::distance(method = dist) %>% stats::hclust(method = hcluter_method)
  clus <- cutree(hc, cuttree)
  d = data.frame(label = names(clus), member = factor(clus))
  map = data.frame(ID = unique(phyloseq::sample_data(ps1_rela)$Group),
                   row.names = unique(phyloseq::sample_data(ps1_rela)$Group),
                   Group = unique(phyloseq::sample_data(ps1_rela)$Group))
  dd = merge(d, map, by = "row.names", all = F)
  row.names(dd) = dd$Row.names
  dd$Row.names = NULL
  p3 = ggtree(hc) %<+% dd + geom_tippoint(size = 5, shape = 21,
                                          aes(fill = member, x = x)) + geom_tiplab(aes(color = member,
                                                                                       x = x * 1.2), hjust = 1)
  p3
  p3 <- p3 + ggnewscale::new_scale_fill()
  head(grotax)
  p4 <- facet_plot(p3, panel = "Stacked Barplot", data = table,
                   geom = ggstance::geom_barh, mapping = aes(x = Abundance,
                                                             fill = !!sym(j)), color = "black", stat = "identity")
  p4
  return(list(p, p1, p3, p4, Taxonomies))
}
