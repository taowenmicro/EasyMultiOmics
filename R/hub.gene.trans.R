#' @title Cluster and Visualize Gene Expression Data
#'
#' @description This function performs clustering on gene expression data and generates visualization plots.
#' It takes a phyloseq object, processes the data, performs clustering using specified method,
#' and returns cluster visualizations and data.
#'
#' @param ps_pro A phyloseq object containing the gene expression data
#' @param id Character string specifying the column name to use for gene IDs (default: "KO_id")
#' @param cluster.method Character string specifying the clustering method to use (default: "mfuzz")
#' @param cluster.num Numeric specifying the number of clusters to generate (default: 6)
#'
#' @return A list containing:
#' \itemize{
#'   \item p1 - Line plot visualization of clusters
#'   \item p2 - Combined heatmap and line plot visualization of clusters
#'   \item dat - Cluster data including cluster assignments
#' }
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njau.edu.cn}
#' @examples
#' \dontrun{
#' # Example usage:
#' result <- hub.gene.trans(ps_pro = my_phyloseq_object)
#' result$p1  # View line plot
#' result$p2  # View combined plot
#' cluster_data <- result$dat  # Access cluster data
#' }
#'
#' @export
hub.gene.trans <- function(ps_pro= ps.trans, id = "KO_id", cluster.method = "mfuzz", cluster.num = 6) {

  # Check taxonomy table columns
  tax = tax_table(ps_pro) %>% head() %>% as.data.frame()

  # Aggregate data by gene ID
  ps_pro1 = ps_pro %>% tax_glom_wt(id) # GeneID

  # Prepare data matrix
  dat <- ps_pro1 %>%
    psmelt() %>%
    dplyr::select(OTU, Sample, Abundance, Group) %>%
    dplyr::group_by(Group, OTU) %>%
    dplyr::summarise(mean = mean(Abundance), .groups = "drop") %>%
    pivot_wider(names_from = Group, values_from = mean) %>% as.data.frame() %>%
    column_to_rownames("OTU")

  # Get clusters
  cm <- clusterData(obj = dat,
                    cluster.method = cluster.method,
                    cluster.num = cluster.num)

  # Create visualizations
  p1 <- visCluster(object = cm,
                   plot.type = "line")

  p2 <- visCluster(object = cm,
                   plot.type = "both",
                   ggplot.panel.arg = c(4, 1, 24, "grey90", NA),
                   ctAnno.col = ggsci::pal_npg()(7))

  # Extract cluster data
  dat = cm$cluster.list

  # Return results
  return(list(p1 = p1, p2 = p2, dat = dat))
}
