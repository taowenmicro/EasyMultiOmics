#' @title Perform reaction pathway enrichment analysis
#'
#' @description
#' This function conducts reaction pathway enrichment analysis based on the specified differential method (wilcox, t.test, etc.) for a phyloseq format file.
#' @param ps A phyloseq format file used as an alternative for the input containing metabolite composition table, metabolite classification table, and sample metadata.
#' @param dif.method The differential method to use for analysis (default is "wilcox").
#' @return
#' A list containing the following components:
#' \describe{
#'   \item{Plots}{List of ggplot objects representing the reaction pathway enrichment analysis plots.}
#'   \item{Plotdata}{List of data frames containing reaction pathway enrichment analysis results.}
#'   }
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njau.edu.cn}
#' @export
#' @examples
#' \dontrun{
#' ps.ms3 = ps.ms %>% tax_glom_wt("KEGGID")
#' res3 = reaction.show.ms(ps = ps.ms3, dif.method = "wilcox")
#' res3$Plots$OE.WT.plot
#' dat = res3$Plotdata$OE.WT
#' head(dat)
#' }

reaction.show.ms = function(
    ps,
    dif.method = "wilcox"
){
  result1 = statSuper(ps = ps,group  = "Group",artGroup = NULL,method = dif.method)
  ps %>% sample_data()
  rank_names(ps)
  group = sample_data(ps)$Group %>% unique()
  com = combn(group,2)
  res = list()
  plot = list()
  # i = 1
  for (i in 1:dim(com)[2]) {
    a = colnames(result1)[str_detect(colnames(result1),"[fdr]")]
    b = colnames(result1)[str_detect(colnames(result1),com[2,i])]
    c = colnames(result1)[str_detect(colnames(result1),com[1,i])]
    d = intersect(a,b)
    e = intersect(d,c)

    tem = e
    n.var= as.name(tem)
    tem2 = result1 %>%
      dplyr::filter(!!n.var < 0.1) %>%
      .$Row.names

    tem2 = tem2[!is.na(tem2)]
    id.tem = tem2 %>% strsplit("[,]") %>%
      sapply( `[`, 1)

    dat = db.compound.reaction.bins
    head(dat)
    colnames(dat)[1] = "compound"
    # 数据库格式调整
    head(dat)
    dat = dat%>%
      separate_rows(reaction, sep = "\\|")

    tax = ps %>%phyloseq::tax_table() %>% as.data.frame()
    head(tax)

    allkeggid <- data.frame(ID = row.names(tax))
    allkeggid$ID = allkeggid$ID%>% strsplit("[,]") %>%
      sapply( `[`, 1)
    head(allkeggid)


    total = dat %>% filter(compound %in% c(paste0("cpd:",allkeggid$ID))) %>% dplyr::select(2,1)
    head(total)


    x <- clusterProfiler::enricher(gene = c(paste0("cpd:",id.tem)),
                                   TERM2GENE = total,minGSSize = 1,pvalueCutoff = 1,qvalueCutoff = 1)

    dat = as.data.frame(x@result)
    head(dat)
    colnames(dat)

    dat = as.data.frame(x@result)
    df <- dat %>%
      arrange(desc(Count)) %>%
      separate(`GeneRatio`,into=c("A","B"),sep="/") %>%
      mutate(A=as.numeric(A),B=as.numeric(B)) %>%
      mutate(count=A/B)  %>% arrange(Count)
    head(df)
    row.names(df) = NULL

    df$Description <- factor(df$Description,levels = c(df$Description %>% as.data.frame() %>% pull()))
    colnames(df)

    df$p.adjust
    df$ID2 = df$ID %>% strsplit("[|]") %>%
      sapply( `[`, 1)

    df$ID2  = gsub("rn:","",df$ID2)
    tab =  db.compound.reaction.dis
    head(tab)

    df2 = df %>% left_join(tab,by = c("ID2" = "V1")) %>%
      distinct(ID2,.keep_all = TRUE) %>%
      arrange(desc(count))
    head(df2)
    res[[i]] = df2
    names(res)[i] = paste(com[2,i],com[1,i],sep = ".")
    tem = df2 %>%
      filter(p.adjust< 0.05)
    head(df2)
    if (dim(tem)[1] == 0) {
      laby = "Top 20 pathways"
      df3 = df2 %>% arrange(desc(p.adjust)) %>%
        head(20)
    } else{
      laby = "enriched pathways"
      df3 = df2 %>%
        filter(p.adjust< 0.05)
    }
    head(df3)

    p = df3 %>%
      ggplot(aes(x = count,y = reorder(V2,count)))+
      geom_bar(aes(x = count,y = reorder(V2,count),fill=pvalue),stat = "identity",width = 0.8)+
      geom_text(aes(x = count,y = reorder(V2,count),label = geneID),hjust = 3)+
      # geom_point(aes(color=pvalue,fill=pvalue),pch=21)+

      scale_color_gradientn(colours = (rev(RColorBrewer::brewer.pal(11,"RdBu"))))+
      scale_fill_gradientn(colours =(rev(RColorBrewer::brewer.pal(11,"RdBu"))))+
      guides(size=guide_legend(title="Count"))+
      labs(x=NULL,y="Pathways")+
      theme(axis.title = element_blank(),
            axis.text.x=element_text(color="black",angle =0,hjust=0.5,vjust=0.5, margin = ggplot2::margin(b =5)),
            axis.text.y=element_text(color="black",angle =0,hjust=1,vjust=0.5),
            panel.background = element_rect(fill = NA,color = NA),
            panel.grid.minor= element_line(size=0.2,color="#e5e5e5"),
            panel.grid.major = element_line(size=0.2,color="#e5e5e5"),
            panel.border = element_rect(fill=NA,color="black",size=1,linetype="solid"),
            legend.key=element_blank(),
            legend.title = element_text(color="black",size=9),
            legend.text = element_text(color="black",size=8),
            legend.spacing.x=unit(0.1,'cm'),
            legend.key.width=unit(0.5,'cm'),
            legend.key.height=unit(0.5,'cm'),
            legend.background=element_blank(),
            legend.box="horizontal",
            legend.box.background = element_rect(color="black"),
            legend.position = c(1,0),legend.justification = c(1,0))+
      scale_y_discrete(labels = function(y) str_wrap(y,width=80))
    plot[[i]] = p
    names(plot)[i] = paste(com[2,i],com[1,i],"plot",sep = ".")

  }
  return(list(plots = plot,plotdata = res))
}
