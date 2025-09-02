#' @title KEGG Pathway Enrichment Analysis for Microbial Communities
#'
#' @description
#' The `KEGG_enrich.metf` function performs KEGG pathway enrichment analysis for microbial communities. It identifies enriched pathways based on differential abundance results and visualizes the results for different group comparisons.
#'
#' @param ps A `phyloseq` object containing OTU (or ASV) table, taxonomic annotations, and sample metadata.
#' @param dif A data frame containing differential abundance results for microbial features (e.g., genes, pathways). Must include columns with group-specific comparisons and significance labels.
#'
#' @return A list of data frames, where each data frame contains KEGG enrichment results for a specific group comparison.
#'
#' @details
#' The function iterates through all pairwise group comparisons, performs KEGG pathway enrichment using the `clusterProfiler::enrichKEGG` function, and returns the enrichment results for each comparison. It filters pathways based on adjusted p-value and q-value thresholds and prepares data for visualization.
#'
#' Steps:
#' \itemize{
#'   \item Extracts all pairwise group comparisons from the sample metadata in the `phyloseq` object.
#'   \item For each comparison, identifies significantly different features (`ko` identifiers).
#'   \item If no significant features are found, selects the top 500 features for analysis.
#'   \item Runs KEGG pathway enrichment analysis using the `enrichKEGG` function.
#'   \item Filters and organizes enrichment results for visualization (e.g., bubble plots).
#' }
#'
#' @examples
#' \dontrun{
#' res = EdgerSuper.metf (ps = ps.kegg,group  = "Group",artGroup = NULL)
#' dat = res[[2]]
#' res2 = KEGG_enrich.metf(ps =  ps.kegg,  diffpath = diffpath,dif = dat)
#'
#' }
#'
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#'
#' @export
KEGG_enrich.metf <- function(ps = ps,
                             # diffpath = diffpath,
                             #  method.diff = md[4],
                             # enrichpath = enrichpath,
                             dif = NULL
){
  getOption("clusterProfiler.download.method")
  R.utils::setOption( "clusterProfiler.download.method",'auto' )
  Desep_group <- ps %>% sample_data() %>%
    .$Group %>%
    as.factor() %>%
    levels() %>%
    as.character()
  Desep_group
  cbtab = combn(Desep_group,2)
  cbtab
  # j= 2
  result_list <- list()
  for (j in 1:dim(cbtab)[2]) {

    Desep_group = cbtab[,j]
    group = paste(Desep_group[1],Desep_group[2],sep = "-")
    print(group)
    #path_ko_plot = paste(enrichpath,"/Ko",group,sep = "")
    #dir.create(path_ko_plot)

    # if (is.null(dif)) {
    #   filename = paste0(diffpath,"/",method.diff,"/","differ.ko.csv")  #paste(diffpath,"/t_all_gene.csv",sep = "")
    #   dif = read_csv(filename)
    # } else if(!is.null(dif)) {
    #   dif = dif
    # }

    dif = dif
    id = paste(group,"level",sep = "")
    id
    head(dif)
    colnames(dif)

    #colnames(dif)[1] = "X1"

    dif$X1 = row.names(dif)

    index <- dif %>%
      dplyr:: select(!!ensym(id),X1) %>%
      dplyr::filter(!!ensym(id) != "nosig") %>%
      .$X1
    num.t = index

    if (length(index) <5) {
      print("nosig")
      index = dif %>%
        select(!!ensym(id),X1) %>%
        # dplyr::filter(!!ensym(id) != "nosig") %>%
        .$X1 %>%
        head(500)

    }

    # ?enrichKEGG
    index=  gsub("ko:","",index)
    kk <- enrichKEGG(gene = index,
                     organism='ko', keyType='kegg',
                     pvalueCutoff = 1)
    head(kk,6)


    plot_GO = subset(kk@result,qvalue < 0.05  & p.adjust < 0.05)

    id.s = nrow(plot_GO)
    if (nrow(plot_GO) < 5) {
      plot_GO  = kk@result
    }
    plot_GO$sig = " "
    plot_GO$sig[plot_GO$pvalue<0.05 &plot_GO$p.adjust <0.2] = "*"

    #-绘制富集气泡图
    ##排序
    plot_GO<- dplyr::arrange(plot_GO,  Count)
    head(plot_GO)


    plot_GO$Description = factor(plot_GO$Description,levels = (plot_GO$Description))

    plot_GO$Description %>% table()

    # p = ggplot(plot_GO,aes( x = Count, y = Description,fill = "#FF0000",size = Count )) +
    #   geom_point(pch = 21) +
    #   geom_text(aes(x = Count, y = Description,label = sig),vjust = 0.6) +
    #   theme_classic()
    # p
    # p
    #
    #
    # p2 = ggplot(plot_GO,aes( x = Count, y = Description,fill = "#FF0000",size = Count )) +
    #   geom_bar(stat = "identity") +
    #   geom_text(aes(x = Count, y = Description,label = sig),vjust = 0.6) +
    #   theme_classic()
    # p2

    result_list[[j]] <- plot_GO
    names( result_list)[j] <-group

  }

  return(result_list)
}
