
#' @title KEGG Enrichment Bubble Plot
#'
#' @description
#' `buplot.micro` creates bubble plots for KEGG enrichment results. It visualizes the relationships between z-scores (effect size), p-values, and the pathway descriptions.
#'
#' @param dt A data frame containing KEGG enrichment results, typically the output from `clusterProfiler::enrichKEGG`. It must include columns for `geneID`, `Description`, `pvalue`, and `p.adjust`.
#' @param id The column name in `dif` corresponding to the group or classification (e.g., "enriched", "depleted").
#'
#' @return A list containing two ggplot2 objects:
#' \itemize{
#'   \item Bubble plot with pathway descriptions labeled.
#'   \item Bubble plot without pathway descriptions labeled.
#' }
#'
#' @details
#' The function calculates a z-score for each pathway based on the proportion of "enriched" and "depleted" genes associated with the pathway, normalized by the pathway size. It then colors the bubbles by the combined z-score and p-value and optionally labels the pathways on the plot.
#'
#' Steps:
#' \itemize{
#'   \item Extract `geneID` from the KEGG result and match it with the `dif` data frame.
#'   \item Calculate the z-score for each pathway as the difference in enriched and depleted counts normalized by the pathway size.
#'   \item Visualize the results using ggplot2.
#' }
#'
#' @examples
#' \dontrun{
#'result = buplot.micro(dt=dat1,id = id)
#'p1 = result[[1]]
#'p2 = result[[2]]
#' }
#'
#' @author Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#'
#' @export
buplot.micro <- function(dt  = kk@result,id = id) {
  data = dt
  data$GeneCount <- as.numeric(sapply(strsplit(data$GeneRatio, "/"), function(x) as.numeric(x[1])))
  data$TotalGenes <- as.numeric(sapply(strsplit(data$GeneRatio, "/"), function(x) as.numeric(x[2])))
  data$HitRatio <- data$GeneCount / data$TotalGenes

  data$level = ifelse(data$pvalue < 0.05, "signif","nosig")

  data$color = data$HitRatio + (-log10(data$pvalue))

  p1 <- ggplot(data, aes(x =HitRatio, y = -log10(pvalue)))  +
    geom_point(pch = 21,aes(fill=color,size =HitRatio))  +
    ggrepel::geom_text_repel(aes(x =HitRatio, y = -log10(pvalue),label= Description)) +
    scale_fill_gradientn(colours =colorRampPalette(c("#F7F4F9","#FFFF33","#FF7F00"))(60)) +
    theme_bw()
  # FileName2 <- paste(path_ko_plot,"/KEGG_enrich_Bubble",group,".pdf", sep = "")
  # ggsave(FileName2, p1, width = 12, height =11,limitsize = FALSE)
  p2 <- ggplot(data, aes(x =HitRatio, y = -log10(pvalue)))  +
    geom_point(pch = 21,aes(fill=color,size =HitRatio ))  +
    # ggrepel::geom_text_repel(aes(x =zscore, y = -log10(pvalue),label= Description)) +
    scale_fill_gradientn(colours =colorRampPalette(c("#F7F4F9","#FFFF33","#FF7F00"))(60)) +
    theme_bw()
  # FileName2 <- paste(path_ko_plot,"/KEGG_enrich_Bubble",group,".pdf", sep = "")
  # ggsave(FileName2, p1, width = 12, height =11,limitsize = FALSE)
  # FileName2 <- paste(path_ko_plot,"/KEGG_enrich_Bubble_nolabel",group,".pdf", sep = "")
  # ggsave(FileName2, p2, width = 8, height =7,limitsize = FALSE)

  return(list(p1,p2))
}

