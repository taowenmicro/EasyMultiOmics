#' hub.network.omics
#'
#' This function calculates the hub scores from a correlation matrix, identifies the top hub nodes,
#' and visualizes the results in a bar plot. It uses the `igraph` package to create a graph from the
#' correlation data and then computes hub scores. The top `tem` hub nodes are extracted and displayed.
#'
#' @param cor A data frame representing a correlation matrix.
#' @param top An integer specifying how many top hub nodes to display (default: `10`).
#'
#' @return A ggplot object representing the bar plot of the top hub nodes.
#'
#' @examples
#' hub_plot <- hub.network.omics(cor = correlation_matrix, top = 10)
#' print(hub_plot)
#'
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}

hub.network.omics <- function(cor, top = 10) {

  # Create an igraph object from the correlation matrix
  igraph <- make_igraph(as.data.frame(cor))

  # Calculate hub scores and select the top nodes
  hub <- hub_score(igraph)$vector %>%
    sort(decreasing = TRUE) %>%
    head(top) %>%
    as.data.frame()

  # Rename the column for hub scores
  colnames(hub) <- "hub_sca"

  # Create a ggplot bar plot for the top hub nodes
  p <- ggplot(hub) +
    geom_bar(aes(x = hub_sca, y = reorder(row.names(hub), hub_sca)), stat = "identity", fill = "#4DAF4A") +
    theme_base()

  # Return the plot
  return(list(p, hub))
}
