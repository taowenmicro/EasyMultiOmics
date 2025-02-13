#' @title KEGG Pathway Enrichment Analysis with Differential Expression
#' @description
#'
#' This function performs KEGG pathway enrichment analysis using KEGG Orthology (KO)
#' terms. It analyzes differential expression results from the `EdgerSuper.trans()` function
#' and identifies significantly enriched KEGG pathways. The function also generates
#' various plots to visualize the KEGG enrichment results, including bubble and bar plots.
#'
#' @param ps A Phyloseq object containing OTU table and sample data. The default is `ps`,
#' which assumes that the Phyloseq object is provided externally.
#' @param method.diff A string specifying the method for differential analysis. The default
#' is `md[4]`, which should be set appropriately based on the method you wish to use.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{plots}: A list of generated plots, including bubble and bar plots for
#'   KEGG pathway enrichment.
#'   \item \code{plotdata}: A list of data used for plotting, including enriched pathway
#'   information.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- KEGG_enrich.trans(ps = ps.trans)
#' }

KEGG_enrich.trans <- function(ps = ps,
                             method.diff = md[4]
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
  res = EdgerSuper.trans (ps = ps,
                          group  = "Group",
                          artGroup = NULL)
  dif = res[[2]]

  # j= 2
  plot1 = list()
  plot2 = list()
  plot3 = list()
  plot4 = list()
  dat = list()

  for (j in 1:dim(cbtab)[2]) {

    Desep_group = cbtab[,j]
    group = paste(Desep_group[1],Desep_group[2],sep = "-")
    print(group)
    id = paste(group,"level",sep = "")
    id
    head(dif)
    dif$X1 = row.names(dif)
    index <- dif %>%
      select(!!ensym(id),X1) %>%
      filter(!!ensym(id) != "nosig") %>%
      .$X1
    num.t = index
    if (length(index) <5) {
      print("nosig")
      index = dif %>%
        select(!!ensym(id),X1) %>%
        # filter(!!ensym(id) != "nosig") %>%
        .$X1 %>%
        head(500)

    }


    index=  gsub("ko:","",index)
    kk <- enrichKEGG(gene = index,
                     organism='ko', keyType='kegg',
                     pvalueCutoff = 1)

    plot_GO = subset(kk@result,qvalue < 0.05  & p.adjust < 0.05)

    id.s = nrow(plot_GO)
    if (nrow(plot_GO) < 5) {
      plot_GO  = kk@result
    }
    plot_GO$sig = " "
    plot_GO$sig[plot_GO$pvalue<0.05 &plot_GO$p.adjust <0.2] = "*"


    plot_GO<- dplyr::arrange(plot_GO,  Count)



    plot_GO$Description = factor(plot_GO$Description,levels = (plot_GO$Description))
    plot_GO$Description %>% table()
    p1 = ggplot(plot_GO,aes( x = Count, y = Description,fill = "#FF0000",size = Count )) +
      geom_point(pch = 21) +
      geom_text(aes(x = Count, y = Description,label = sig),vjust = 0.6) +
      theme_classic()
    p1

    plot1[[j]] = p1
    names(plot1)[j] = group

    # if (id.s < 5) {
    #   FileName2 <- paste(path_ko_plot,"/KEGG_all_point.nosig",group,".pdf", sep = "")
    #   # , device = cairo_pdf, family = "Times New Roman"
    #   ggsave(FileName2, p, width = 18, height =2* (nrow(plot_GO)/3),limitsize = FALSE)
    #   FileName2 <- paste(path_ko_plot ,"/KEGG_all_point.nosig",group,".csv", sep = "")
    #   # , device = cairo_pdf, family = "Times New Roman"
    #   write_csv(plot_GO,FileName2)
    # } else{
    #   FileName2 <- paste(path_ko_plot ,"/kegg_all_point.",group,".pdf", sep = "")
    #   # , device = cairo_pdf, family = "Times New Roman"
    #   ggsave(FileName2, p, width = 18, height =2* (nrow(plot_GO)/3),limitsize = FALSE)
    #   FileName2 <- paste(path_ko_plot,"/kegg_all_point.",group,".csv", sep = "")
    #   # , device = cairo_pdf, family = "Times New Roman"
    #   write_csv(plot_GO,FileName2)
    # }


    p2 = ggplot(plot_GO,aes( x = Count, y = Description,fill = "#FF0000",size = Count )) +
      geom_bar(stat = "identity") +
      geom_text(aes(x = Count, y = Description,label = sig),vjust = 0.6) +
      theme_classic()
    p2

    plot2[[j]] = p2
    names(plot2)[j] = group

    dat[[j]] = plot_GO
    names(dat)[j] = group
    # if (id.s < 5) {
    #
    #   FileName2 <- paste(path_ko_plot,"/KEGG_all_bar.nosig",group,".pdf", sep = "")
    #
    #   # , device = cairo_pdf, family = "Times New Roman"
    #   ggsave(FileName2, p, width = 18, height =2* (nrow(plot_GO)/3),limitsize = FALSE)
    #
    #   FileName2 <- paste(path_ko_plot ,"/KEGG_all_bar.nosig",group,".csv", sep = "")
    #   # , device = cairo_pdf, family = "Times New Roman"
    #   write_csv(plot_GO,FileName2)
    #
    # } else{
    #   FileName2 <- paste(path_ko_plot,"/kegg_all_bar.",group,".pdf", sep = "")
    #   # , device = cairo_pdf, family = "Times New Roman"
    #   ggsave(FileName2, p, width = 18, height =2* (nrow(plot_GO)/3),limitsize = FALSE)
    #   FileName2 <- paste(path_ko_plot,"/kegg_all_bar.",group,".csv", sep = "")
    #   # , device = cairo_pdf, family = "Times New Roman"
    #   write_csv(plot_GO,FileName2)
    # }

    if (length(num.t) < 5) {

    } else{

      #-绘制富集气泡图
      print("1")
      ?buplot
      result = buplot.micro(dt  = plot_GO,id = id)
      print("1")

      p3 = result[[1]]
      p4 = result[[2]]


      plot3[[j]] = p3
      names(plot3)[j] = group

      plot4[[j]] = p4
      names(plot4)[j] = group
    }

  }
  return(list(plots = list(plot1,plot2,plot3,plot4),plotdata = dat))

}
