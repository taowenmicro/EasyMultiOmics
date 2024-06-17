



#-------PCA载荷挑选#---------
# pcapath = paste(repath,"/loadingPCA/",sep = "")
# dir.create(pcapath)
# res = loadingPCA(ps = ps)
#
# p = res[[1]]
# p
# dat = res[[2]]
#
# filemane = paste(pcapath,"/PCALoading.pdf",sep = "")
# ggsave(filemane, p, width = 8, height = 6)
#
# FileName <- paste(pcapath,"/Loadsing_pca.csv", sep = "")
# write.csv(dat,FileName,sep = "")

loadingPCA.micro = function(ps = ps,Top = 20){
  count =ps %>%
    ggClusterNet::vegan_otu() %>% t()
  count[is.na(count)] = 0
  norm = t(t(count)/colSums(count,na=TRUE))# * 100 # normalization to total 100
  otu.pca <- stats::prcomp(t(norm), scale= FALSE)

  #提取
  yangpin<-otu.pca$x
  yangpin=as.data.frame(yangpin)
  yangpin$SampleType= phyloseq::sample_data(ps)$Group
  #提取荷载坐标
  bianliang<-otu.pca$rotation
  bianliang=as.data.frame(bianliang)
  head(bianliang)
  dim(norm)
  index = merge(norm ,bianliang, by="row.names",all=F)
  head(index)
  row.names(index)=index$Row.names
  index$Row.names=NULL
  head(index)
  index$id = row.names(index)
  ##手动选择10个最终要的变量 PCA载荷矩阵挑选37个成分提取差异.txt
  index$PCone = index$PC1^2
  top = index %>% arrange(desc(PCone)) %>%
    head(Top)

  head(top)
  top$ID  = top$id
  head(top)

  p=ggplot(top, aes(x = PCone, y = reorder(ID,PCone)))  +
    geom_segment(aes(yend=ID),xend=0,size=3,colour = "#1B9E77" )+
    geom_point(size=4,pch=20, colour = "#1B9E77")+theme_bw()+
    theme(axis.text.x = element_text(colour = "black",size = 20,face = "bold"),
          axis.text.y = element_text(colour = "black",size = 10,face = "bold"))
  p
  return(list(p,index))
}


