
GSVA.metf= function(ps = ps,
                     #  GSVApath = GSVApath,
                      FC = FALSE,
                      adj = FALSE,
                      lg2FC = 0.3,
                      padj = 0.05
){

  kegg <- clusterProfiler::download_KEGG('ko')
  PATH2ID <- kegg $KEGGPATHID2EXTID
  PATH2NAME <- kegg$KEGGPATHID2NAME
  head(PATH2NAME)
  PATH_ID_NAME <- merge(PATH2ID, PATH2NAME, by='from')
  colnames(PATH_ID_NAME) <- c('KEGGID', 'KO', 'DESCRPTION')
  head(PATH_ID_NAME )
  ##查看有多少个ko号
  length(unique(PATH_ID_NAME$KEGGID))
  #转化成list
  group <- split(PATH_ID_NAME$KO,PATH_ID_NAME$KEGGID)
  a <- PATH_ID_NAME$KO %>% unique()
  tab = ps %>% vegan_otu() %>%
    t()
  map = ps %>% sample_data()
  Desep_group <- as.character(levels(as.factor(map$Group)))
  aaa = combn(Desep_group,2)
  row.names(tab) = gsub("ko:","",row.names(tab))
  str(tab)

  topMatrixGSVA <- GSVA::gsva(as.matrix(tab), group,
                        min.sz=10, max.sz=999999,
                        abs.ranking=FALSE,
                        verbose=TRUE)

  topMatrixGSVA = topMatrixGSVA %>% as.data.frame()
  topMatrixGSVA$ID = row.names(topMatrixGSVA)


  tax = PATH_ID_NAME %>%
    dplyr::distinct(KEGGID,.keep_all = TRUE) %>%
    dplyr::inner_join(topMatrixGSVA,by = c("KEGGID" = "ID")
    )

  tax$ID = paste(tax$KEGGID,tax$DESCRPTION,sep = ":")
  row.names(tax) = tax$ID
  dim(tax)
  map$ID = row.names(map)
  heat <- t(scale(t(tax[map$ID])))
  myCol <- colorRampPalette(c("dodgerblue", "black", "yellow"))(100)
  myBreaks <- seq(-1.5, 1.5, length.out=101)
  # filename3 = paste(GSVApath,"/GSVE_all",".csv",sep = "")
  # write.csv(tax,filename3)
  head(heat)

  result <- heatmap_GSVA(heat = heat,
                         map = map,
                         col_cluster =  TRUE,
                         row_cluster =  TRUE,
                         label =  TRUE)
  p1 = result[[1]]


  p2 = result[[2]]

  p1

  plot1 = list( heat, p1, p2)
  # if (dim(heat)[2] < 10) {
  #   w = 2
  # } else{
  #   w = 0.5
  # }
  #
  # filename3 = paste(GSVApath,"/GSVE_all_heatmap",".pdf",sep = "")
  # ggsave(filename3,p1,width = dim( heat)[2]*w, height = dim( heat)[1]/8,limitsize = FALSE)
  # filename3 = paste(GSVApath,"/GSVE_all_bubb",".pdf",sep = "")
  # ggsave(filename3,p2,width = dim( heat)[2]*w, height = dim( heat)[1]/8,limitsize = FALSE)

  design =model.matrix(~ 0 + map$Group)
  colnames(design)=levels(as.factor(map$Group))
  fit <- limma::lmFit(topMatrixGSVA[,map$ID], design)
  fit <- limma::eBayes(fit)
  # i = 1
  plot_list = list()
  plotlist2 = list()
  for (i in 1:dim(aaa)[2]) {

    sds = paste(aaa[1,i],aaa[2,i],sep = "-")
    print("1")
    aa = function(sds = sds,design = design){
      contrast.matrix<- limma::makeContrasts(contrasts =  sds,levels= design)
    }

    contrast.matrix = aa(sds = sds,design = design)


    # contrast.matrix = matrix(nrow = c(2),ncol = 1,c(1,-1))
    colnames( contrast.matrix) = sds

    print(i)
    id = sds

    # contrast.matrix <- divgro(id = sds)
    fit2 <- limma::contrasts.fit(fit, contrast.matrix) %>%
      limma::eBayes()
    #"fdr"(equivalent to "BH");lfc(log2-fold-change )
    results<- limma::decideTests(fit2,
                                 method="global",
                                 adjust.method="BH",
                                 p.value=0.05,
                                 lfc=0)
    summary(results)
    print("1")
    x<- limma::topTable(fit2,coef=1, number=10000, adjust.method="BH", sort.by="p",resort.by=NULL)
    x <- data.frame(rownames(x), x)
    colnames(x) <- c("Pathway","Log2FoldChange","MeanExpression","tStat","Pvalue","AdjustedPval","Bvalue")
    # head(x)
    # padj  = 0.05

    if (adj == TRUE) {
      WT <-subset(x,AdjustedPval < padj )
      head(WT)
    } else {
      WT <-subset(x,Pvalue < padj )
      head(WT)
    }
    # lg2FC = 0.3
    if (FC == TRUE) {
      WT <-subset(WT,abs(Log2FoldChange) > lg2FC )
      head(WT)
    } else {
      WT = WT
    }
    dim(WT)
    WT = WT %>% dplyr::arrange(Log2FoldChange) %>% as.data.frame()
    row.names(WT) = WT$Pathway
    #Filter the GSVA object to only include significant pathways
    head(tax)
    sigtab <- tax %>%
      dplyr::filter(KEGGID  %in% rownames(WT))
    row.names(sigtab) = sigtab$ID
    sigtab <- sigtab[match(rownames(WT),sigtab$KEGGID),]
    #Set colour for heatmap
    # require(RColorBrewer)
    myCol <- colorRampPalette(c("dodgerblue", "black", "yellow"))(100)
    myBreaks <- seq(-1.5, 1.5, length.out=101)
    heat <- t(scale(t(sigtab[,map$ID]))) %>% as.data.frame()
    heat$id = row.names(heat)
    # sub_heat<- dplyr::filter(heat,id %in% row.names(WT))
    sub_map<- dplyr::filter(as.tibble(map),Group %in% aaa[,i])
    # sub_map$ID = sub_map$
    vars<- c(as.character(sub_map$ID),"id")
    sub_heat = dplyr::select(heat,one_of(vars))
    dim(sub_heat)
    row.names(sub_heat) = sub_heat$id
    sub_heat$id = NULL
    sub_heat = as.matrix(sub_heat)

    # library(gplots)
    #
    # filename3 = paste(GSVApath,"/GSVE_",paste(aaa[1,i],aaa[2,i],sep = "-"),".csv",sep = "")
    # write.csv(x,filename3)
    #
    head(sub_heat)
    plot_list[[i]] = sub_heat
    names(plot_list)[[i]] = sds


    sub_heat = sub_heat %>% as.data.frame() %>% rownames_to_column("ID") %>%
    #   dplyr::filter(!str_detect(ID,"ko00998")) %>%
      column_to_rownames("ID") %>% as.matrix()


    if (dim(sub_heat)[1]> 2) {

      result <- heatmap_GSVA(heat = sub_heat,
                             map = sub_map,
                             col_cluster =  F,
                             row_cluster =  F,
                             label =  TRUE)
      p1 = result[[1]]

      plotlist2[[i]] = p1
      names(plotlist2)[[i]] = sds


      # dev.off()
      p1
      #
      # if (dim(heat)[2] < 10) {
      #   w = 1.5
      # } else{
      #   w = 1
      # }
      #
      # filename3 = paste(GSVApath,"/GSVE_all_heatmap",paste(aaa[1,i],aaa[2,i],sep = "-"),".pdf",sep = "")
      # ggsave(filename3,p1,width = dim( sub_heat)[2]*w, height = dim( sub_heat)[1]/6,limitsize = FALSE)
      # filename3 = paste(GSVApath,"/GSVE_all_bubb",paste(aaa[1,i],aaa[2,i],sep = "-"),".pdf",sep = "")
      # ggsave(filename3,p2,width = dim( sub_heat)[2]*w, height = dim(sub_heat)[1]/6,limitsize = FALSE)
      #
    }

  }
  return(list(plot1, plot_list,plotlist2))
}
