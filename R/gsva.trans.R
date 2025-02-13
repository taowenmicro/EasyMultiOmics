# res = gsva.trans(ps = ps.trans,FC = FALSE,
#                  adj = FALSE,lg2FC = 0.3,
#                  padj = 0.05)
# #
# res$plots[[1]][[1]]
# #
# res$plots[[2]]$`KO OE`



gsva.trans  = function(ps = ps,
                       FC = FALSE,
                       adj = FALSE,
                       lg2FC = 0.3,
                       padj = 0.05
){
  getOption("clusterProfiler.download.method")
  R.utils::setOption( "clusterProfiler.download.method",'auto' )

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
  head(tab)
  map = ps %>% sample_data()
  Desep_group <- as.character(levels(as.factor(map$Group)))
  aaa = combn(Desep_group,2)

  row.names(tab) = gsub("ko","",row.names(tab))
  topMatrixGSVA <- gsva(as.matrix(tab),group, min.sz=10, max.sz=999999, abs.ranking=FALSE, verbose=TRUE)

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
  dat = tax
  result <- heatmap_GSVA(heat = heat,
                         map = map,
                         col_cluster =  TRUE,
                         row_cluster =  TRUE,
                         label =  TRUE)
  p1 = result[[1]]
  p2 = result[[2]]
  plot1 = list()
  plot1[[1]] = p1
  plot1[[2]] = p2

  # if (dim(heat)[2] < 10) {
  #   w = 2
  # } else{
  #   w = 0.5
  # }
  # filename3 = paste(GSVApath,"/GSVE_all_heatmap",".pdf",sep = "")
  # ggsave(filename3,p1,width = dim( heat)[2]*w, height = dim( heat)[1]/8,limitsize = FALSE)
  # filename3 = paste(GSVApath,"/GSVE_all_bubb",".pdf",sep = "")
  # ggsave(filename3,p2,width = dim( heat)[2]*w, height = dim( heat)[1]/8,limitsize = FALSE)


  design =model.matrix(~ 0 + map$Group)
  colnames(design)=levels(as.factor(map$Group))
  fit <- limma::lmFit(topMatrixGSVA[,map$ID], design)
  fit <- limma::eBayes(fit)
  # i = 1

  plot2 = list()
  plot3 = list()
  plot4 = list()
  dat2 = list()

  aa = function(sds = sds,design = design){
    contrast.matrix<- limma::makeContrasts( contrasts = sds,levels= design)
  }
  contrast.matrix = aa(sds = sds,design = design)

  for (i in 1:dim(aaa)[2]) {

    sds = paste(aaa[1,i],aaa[2,i],sep = "-")
    print("1")


    # contrast.matrix = matrix(nrow = c(2),ncol = 1,c(1,-1))
    # colnames( contrast.matrix) = sds

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

    if (adj) {
      WT <-subset(x,AdjustedPval < padj )
      head(WT)
    } else {
      WT <-subset(x,Pvalue < padj )
      head(WT)
    }
    # lg2FC = 0.3
    if (FC) {
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
    sub_map <- dplyr::filter(as.tibble(map),Group %in% aaa[,i])
    # sub_map$ID = sub_map$
    vars<- c(as.character(sub_map$ID),"id")
    sub_heat = dplyr::select(heat,one_of(vars))
    head(sub_heat)



    row.names(sub_heat) = sub_heat$id
    sub_heat$id = NULL
    sub_heat = as.matrix(sub_heat)



    otu_table = sub_heat
    head(otu_table)

    design = as.data.frame(sample_data(ps))
    ## 计算相对丰度，计算每个物种丰度均值，按照均值排序
    OTU = as.matrix(otu_table)
    norm = t(t(OTU)/colSums(OTU,na=TRUE)) #* 100 # normalization to total 100
    norma = norm %>%
      t() %>% as.data.frame()
    #数据分组计算平均值
    iris.split <- split(norma,as.factor(sub_map$Group))

    iris.apply <- lapply(iris.split,function(x)colMeans(x,na.rm = TRUE))
    # 组合结果
    norm2 <- do.call(rbind,iris.apply)%>%
      t()
    norm2 = as.data.frame(norm2)
    norm2$mean=apply(norm2,1,mean)
    norm2$ID = row.names(norm2)
    colnames(norm2)
    abun = norm2
    head(abun)
    abun$mean = NULL

    abun_a = reshape2::melt(abun,
                            id.var = c("ID"),
                            variable.name = "id",
                            value.name = "count")


    head(abun_a)

    abun_a$iid = rep(paste(1:(length(abun_a$id)/2)),2)
    id.all = abun_a$iid %>% unique()
    for (ii in id.all) {
      tem = abun_a %>%
        dplyr::filter(iid == ii) %>% .$count
      if (tem[1]*tem[2] >= 0) {
        abun_a$count[abun_a$iid == ii] = abs(abun_a$count[abun_a$iid == ii])
      }

      if (tem[1]*tem[2] < 0) {
        abun_a$count[abun_a$iid == ii] = tem + min(tem) %>% abs()
      }
    }


    library(plyr)
    abun_a1 = ddply(abun_a,"iid",transform,percount = count/sum(count)*100)
    head(abun_a1)
    mi = c( "firebrick3","navy")

    head(abun_a1)
    tem2 = abun_a1 %>%
      distinct(iid, .keep_all = TRUE) %>%
      arrange(percount) %>%.$ID

    abun_a1$ID = factor(abun_a1$ID,levels = as.character(tem2))
    p3 = abun_a1  %>%
      ggplot(aes(x = ID, y = percount, fill =id, group =id)) +
      geom_bar(stat = 'identity') + scale_fill_manual(values = mi)+
      # theme_classic()+
      theme_classic()+
      theme(axis.text.x=element_text(angle=90,vjust=0.5, hjust=1,size = 10))
    p3
    plot2[[i]] = p3
    names(plot2)[i] = paste(aaa[1,i],aaa[2,i])
    # filename3 = paste(GSVApath,"/GSVE_all_Abuncance.relative",paste(aaa[1,i],aaa[2,i],sep = "-"),".pdf",sep = "")
    # ggsave(filename3,p,width = dim( sub_heat)[2]/4, height = 6,limitsize = FALSE)

    sub_heat = sub_heat %>% as.data.frame() %>% rownames_to_column("ID") %>%
      filter(!str_detect(ID,"ko00998")) %>%
      column_to_rownames("ID") %>% as.matrix()
    if (dim(sub_heat)[1]> 2) {

      result <- heatmap_GSVA(heat = sub_heat,
                             map = sub_map,
                             col_cluster =  F,
                             row_cluster =  F,
                             label =  TRUE)
      p4 = result[[1]]
      p5 = result[[2]]
      # dev.off()
      plot3[[i]] = p4
      names(plot3)[i] = paste(aaa[1,i],aaa[2,i])

      plot4[[i]] = p5
      names(plot4)[i] = paste(aaa[1,i],aaa[2,i])
      dat2[[i]] = sub_heat
      names(dat2)[i] = paste(aaa[1,i],aaa[2,i])
    }

  }
  return(list(plots = list(plot1,plot2,plot3,plot4),plotdata = dat,plotdsata.g = dat2))
}
