
#' @title Multi-Group Volcano Plot for Differential Features
#'
#' @description
#' The `Mui.Group.volcano.metm` function generates a multi-group volcano plot
#' to visualize the differential features (e.g., taxa or genes) across multiple groups.
#' It highlights significantly enriched or depleted features for each group
#' and displays the top-ranked features based on their fold changes.
#'
#' @param res A data frame containing the differential analysis results.
#' Must include the following columns:
#' \itemize{
#'   \item `logFC`: Log2 fold change values.
#'   \item `level`: Significance levels (e.g., "enriched", "depleted", or "nosig").
#'   \item `group`: Group or cluster identifier.
#' }
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item  A fully styled multi-group volcano plot with fold changes, top-ranked features, and group annotations.
#'   \item  A basic multi-group volcano plot with group annotations.
#'   \item  The input data frame with an additional `ID` column for unique identifiers.
#'   \item  A data frame containing the top 5 features (based on absolute fold changes) for each group.
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Computes the top 5 features with the highest absolute fold changes for each group.
#'   \item Creates background columns for each group based on the maximum and minimum fold changes.
#'   \item Adds jittered scatter points to represent individual features' fold changes.
#'   \item Annotates the top-ranked features for each group.
#'   \item Styles the volcano plot with group-specific color tiles and customized aesthetics.
#' }
#'
#' @examples
#' \dontrun{
#' res =  EdgerSuper2.micro (ps = ps.16s,group  = "Group",artGroup =NULL, j = "OTU")
#' res2 = Mui.Group.volcano.metm(res = res)
#' }
#'
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @export


Mui.Group.volcano.metm =function (res = res)
{
  res$ID = row.names(res)
  datv = res
  tm.g <- function(data) {
    id = data$group %>% unique()
    for (i in 1:length(id)) {
      tem = filter(data, group == id[i], level != "nosig") %>%
        distinct(ID, .keep_all = TRUE) %>% top_n(5, abs(logFC))
      if (i == 1) {
        tem2 = tem
      }
      else {
        tem2 = rbind(tem2, tem)
      }
    }
    return(tem2)
  }
  top <- tm.g(datv)
  head(datv)
  tem = datv %>% group_by(group) %>% dplyr::summarise(max = max(logFC),
                                                      min = min(logFC)) %>% as.data.frame()
  col1 <- data.frame(x = tem$group, y = tem$max)
  col2 <- data.frame(x = tem$group, y = tem$min)
  p1 <- ggplot() + geom_col(data = col1, mapping = aes(x = x,
                                                       y = y), fill = "#dcdcdc", alpha = 0.6) + geom_col(data = col2,
                                                                                                         mapping = aes(x = x, y = y), fill = "#dcdcdc", alpha = 0.6)
  p1
  head(datv)
  p2 <- ggplot() + geom_col(data = col1, mapping = aes(x = x,
                                                       y = y), fill = "#dcdcdc", alpha = 0.4) + geom_col(data = col2,
                                                                                                         mapping = aes(x = x, y = y), fill = "#dcdcdc", alpha = 0.4) +
    geom_jitter(data = datv, aes(x = group, y = logFC, color = level),
                size = 1, width = 0.4) + scale_color_manual(name = NULL,
                                                            values = c("#4393C3", "#FC4E2A", "grey40")) + labs(x = "",
                                                                                                               y = "log2(FoldChange)")
  p2
  dfcol <- data.frame(x = tem$group, y = 0, label = tem$group)
  dfcol$group <- tem$group
  library(RColorBrewer)
  library(MetBrewer)
  tile_color <- met.brewer("Thomas", length(tem$group))
  p3 <- p2 + geom_tile(data = dfcol, aes(x = x, y = y), height = 1.75,
                       color = "black", fill = tile_color, alpha = 0.6, show.legend = F) +
    ggrepel::geom_text_repel(data = dfcol, aes(x = x, y = y,
                                               label = group), size = 3.5, color = "white") + theme_classic()
  p3
  p4 <- p3 + geom_text_repel(data = top, aes(x = group, y = logFC,
                                             label = ID), force = 1.2, arrow = arrow(length = unit(0.008,
                                                                                                   "npc"), type = "open", ends = "last"))
  p4
  p5 <- p4 + theme_minimal() + theme(axis.title = element_text(size = 18,
                                                               color = "black", face = "bold"), axis.line.y = element_line(color = "black",
                                                                                                                           size = 1.2), axis.line.x = element_blank(), axis.text.x = element_blank(),
                                     panel.grid = element_blank(), legend.position = "top",
                                     legend.direction = "vertical", legend.justification = c(1,
                                                                                             0), legend.text = element_text(size = 12))
  p5
  return(list(p5, p3, datv, top))
}
