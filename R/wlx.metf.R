wlx.metf = function(
    ps = NULL,
    artGroup = NULL,
    group  = "Group",
    pvalue = 0.05
){

  count  = ps %>%
    transform_sample_counts( function(x) x / sum(x) * 10000 ) %>% vegan_otu() %>%
    t()
  Desep_group <- ps %>% sample_data() %>%
    .$Group %>%
    as.factor() %>%
    levels() %>%
    as.character()
  if ( is.null(artGroup)) {
    #--构造两两组合#-----
    aaa = combn(Desep_group,2)
    # sub_design <- as.data.frame(sample_data(ps))
  } else if (!is.null(artGroup)) {
    aaa  = as.matrix(b )
  }

  plot_list <- list()

  for (j in 1:dim(aaa)[2]) {
    Desep_group1 = aaa[,j]

    group = paste(Desep_group1[1],Desep_group1[2],sep = "-")
    group
    sub_design <- ps %>% sample_data()
    sub_design$ID = row.names(sub_design)
    df_filter<- dplyr::filter(as.tibble(sub_design ),Group %in% Desep_group1)
    head(df_filter)
    a = as.data.frame(count)
    a = a[as.character(df_filter$ID)]
    rep = length(as.character(df_filter$ID))/2
    a = as.matrix(a)

    x <- wilcox(a = a,rep = rep)
    head(x)
    x$level = as.factor(ifelse(as.vector(x$fdr) < 0.05 & x$log2_FC > 0 , "enriched",
                               ifelse(as.vector(x$fdr) < 0.05 &x$log2_FC < 0, "depleted","nosig")))

    x = data.frame(row.names = row.names(a),logFC = as.numeric(x$log2_FC),
                   level = x$level,p = as.numeric(x$Pvalue))
    head(x)
    x$logFC
    p <- ggplot(x,aes(x =logFC ,y = -log2(p), colour=level)) +
      geom_point() +
      geom_hline(yintercept=-log10(0.2),
                 linetype=4,
                 color = 'black',
                 size = 0.5) +
      geom_vline(xintercept=c(-1,1),
                 linetype=3,
                 color = 'black',
                 size = 0.5) +
      # geom_text_repel(data=x, aes(x=logFC, y=-log2(p_adjust), label=row.names(x)), size=1) +
      scale_color_manual(values = c('blue2','red2', 'gray30')) +
      ggtitle(group) + theme_bw()

    p

    # file = paste(path,"/",group,"t.test_Volcano.pdf",sep = "")
    # ggsave(file,p,width = 8,height = 6)
    colnames(x) = paste(group,colnames(x),sep = "")
    plot_list[[paste(group, "plot", sep = "_")]] <- p


    if (j == 1) {
      table =x
    } else if (j != 1) {
      table = cbind(table,x)
    }

  }
  #data_list[["combined_table"]] <- table
  x = table
  head(x)
  count = as.matrix(count)
  norm = t(t(count)/colSums(count)) #* 100 # normalization to total 100
  dim(norm)
  norm1 = norm %>%
    t() %>% as.data.frame()
  # head(norm1)
  #数据分组计算平均值
  # library("tidyverse")
  head(norm1)
  iris.split <- split(norm1,as.factor(sub_design$Group))
  iris.apply <- lapply(iris.split,function(x)colMeans(x))
  iris.combine <- do.call(rbind,iris.apply)
  norm2= t(iris.combine) %>%
    as.data.frame()
  x = cbind(x,norm2)
  return(list(plot_list,x))
}
wilcox <- function(a = a,rep = rep){
  Pvalue<-c(rep(0,nrow(a)))
  fdr<-c(rep(0,nrow(a)))
  log2_FC<-c(rep(0,nrow(a)))

  for(i in 1:nrow(a)){
    if(sd(a[i,(1:rep)])==0&&sd(a[i,(rep+1):(rep*2)])==0){
      Pvalue[i] <-"NA"
      log2_FC[i]<-"NA"
    }else{
      y=wilcox.test(as.numeric(a[i,(1:rep)]),as.numeric(a[i,(rep+1):(rep*2)]),exact=FALSE)
      Pvalue[i]<-y$p.value
      log2_FC[i]<-log2((mean(as.numeric(a[i,(1:rep)]))+0.001)/(mean(as.numeric(a[i,(rep+1):(rep*2)]))+0.001))
      fdr[i]=p.adjust(Pvalue[i], "BH")
    }
  }

  out<-cbind(log2_FC,Pvalue,fdr) %>% as.data.frame()
  return(out)

}
