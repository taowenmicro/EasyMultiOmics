#' @title Ordination analysis and plotting
#' @description Input metabolite table, metadata and tree or phyloseq object;
#' support 47 distance type (bray, unifrac, wunifrac ...),
#' 8 ordination method (PCoA, NMDS, ...);
#' output ggplot2 figure, data and statistical test result.
#' @param otu Metabolite composition table.
#' @param map Sample metadata.
#' @param dist Distance type, including "unifrac","wunifrac","manhattan","euclidean","bray","jaccard"... .
#' @param group Column name for groupID in map table(sample metadata).
#' @param method DCA, CCA, RDA, NMDS, MDS, PCoA, PCA, LDA, t-sne.
#' @param pvalue.cutoff Pvalue threshold.
#' @param tax  Metabolite classification table.
#' @param ps A phyloseq format file used as an alternative for the input containing metabolite, tax, and sample metadata.
#' @param pair A logical value for whether to perform analysis between any two groups.
#' @param Micromet statistics by adonis/anosim/MRPP;
#' @return A list containing the following components:
#' \item{p2}{Standard plot of the ordination results.}
#' \item{points}{Data frame with coordinates of samples in the ordination space.}
#' \item{p3}{Annotated plot with sample labels.}
#' \item{pairResult}{Results of pairwise tests if applicable.}
#' \item{title1}{Overall inspection results by anosim or adonis.}
#' \item{eig}{Eigenvalues or variance explained by each axis.}
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' result = ordinate.ms(ps = ps.ms, group = "Group", dist = "bray",
#' method = "PCoA", Micromet = "anosim", pvalue.cutoff = 0.05,pair = F)
#' p3_1 = result[[1]]
#' p3_1
#' p3_2 = result[[3]]
#' p3_2
#' plotdata =result[[2]]
#' head(plotdata)
#' @export

ordinate.ms = function(otu = NULL,tax = NULL,map = NULL,ps = NULL,
         group = "Group", dist = "bray", method ="PCoA",
         Micromet = "adonis", pvalue.cutoff = 0.05,pair=TRUE){
  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)
  # 求取相对丰度#----
  ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) )


  # 排序方法选择#----

  #---------如果选用DCA排序
  if (method == "DCA") {
    # method = "DCA"
    ordi = phyloseq::ordinate(ps1_rela, method=method, distance=dist)
    #提取样本坐标
    points = ordi$rproj[,1:2]
    colnames(points) = c("x", "y") #命名行名
    #提取特征值
    eig = ordi$evals^2
  }

  #---------CCA排序#----
  if (method == "CCA") {
    # method = "CCA"
    ordi = phyloseq::ordinate(ps1_rela, method=method, distance=dist)
    #样本坐标,这里可选u或者v矩阵
    points = ordi$CA$v[,1:2]
    colnames(points) = c("x", "y") #命名行名
    #提取特征值
    eig = ordi$CA$eig^2
  }

  #---------RDA排序#----
  if (method == "RDA") {
    # method ="RDA"
    ordi = phyloseq::ordinate(ps1_rela, method=method, distance=dist)
    #样本坐标,这里可选u或者v矩阵
    points = ordi$CA$v[,1:2]
    colnames(points) = c("x", "y") #命名行名
    #提取特征值
    eig = ordi$CA$eig
  }

  #---------DPCoA排序#----

  #---------MDS#----
  if (method == "MDS") {
    # method = "MDS"
    # ordi = ordinate(ps1_rela, method=ord_meths[i], distance=dist)
    ordi = phyloseq::ordinate(ps1_rela, method=method, distance=dist)
    #
    points = ordi$vectors[,1:2]
    colnames(points) = c("x", "y")
    #
    eig = ordi$values[,1]
  }

  #---------PCoA#----
  if (method == "PCoA") {
    # method = "PCoA"
    unif = phyloseq::distance(ps1_rela , method=dist, type="samples")
    pcoa = stats::cmdscale(unif, k=2, eig=TRUE) # k is dimension, 3 is recommended; eig is eigenvalues
    points = as.data.frame(pcoa$points) # 获得坐标点get coordinate string, format to dataframme
    colnames(points) = c("x", "y") #命名行名
    eig = pcoa$eig
  }

  #----PCA#----

  otu_table = as.data.frame(t(ggClusterNet::vegan_otu(ps1_rela )))
  # head(otu_table)
  # method = "PCA"
  if (method == "PCA") {
    otu.pca = stats::prcomp(t(otu_table), scale.default = TRUE)
    #
    points = otu.pca$x[,1:2]
    colnames(points) = c("x", "y") #命名行名

    eig=otu.pca$sdev
    eig=eig*eig
  }

  # method = "LDA"
  #---------------LDA#----
  if (method == "LDA") {
    # library(MASS)
    data = t(otu_table)
    # head(data)
    data = as.data.frame(data)
    # data$ID = row.names(data)
    data = scale(data, center = TRUE, scale = TRUE)
    dim(data)
    data1 = data[,1:10]
    map = as.data.frame(sample_data(ps1_rela))
    model = MASS::lda(data, map$Group)

    ord_in = model
    axes = c(1:2)
    points = data.frame(predict(ord_in)$x[, axes])
    colnames(points) = c("x", "y") #
    #
    eig= ord_in$svd^2
  }

  #---------------NMDS#----
  # method = "NMDS"

  if (method == "NMDS") {
    #---------
    # i = 5
    # dist = "bray"
    ordi = phyloseq::ordinate(ps1_rela, method=method, distance=dist)
    #样本坐标,
    points = ordi$points[,1:2]
    colnames(points) = c("x", "y") #
    #提取stress
    stress = ordi$stress
    stress= paste("stress",":",round(stress,2),sep = "")
  }


  #---------------t-sne#----
  # method = "t-sne"
  if (method == "t-sne") {
    data = t(otu_table)
    # head(data)
    data = as.data.frame(data)
    # data$ID = row.names(data)
    #
    data = scale(data, center = TRUE, scale = TRUE)

    dim(data)
    map = as.data.frame(sample_data(ps1_rela))
    row.names(map)
    #---------tsne
    # install.packages("Rtsne")
    # library(Rtsne)

    tsne = Rtsne::Rtsne(data,perplexity = 3)

    #
    points = as.data.frame(tsne$Y)
    row.names(points) =  row.names(map)
    colnames(points) = c("x", "y") #命名行名
    stress= NULL
  }



  g = sample_data(ps)$Group %>% unique() %>% length()
  n = sample_data(ps)$Group%>% length()
  o = n/g

  if (o >= 3) {
    title1 =MicroTest.ms(ps = ps1_rela, Micromet = Micromet, dist = dist)
    title1
  } else{
    title1 = NULL
  }



  if (pair==TRUE) {
    #----两两比较#----
    pairResult = pairMicroTest.ms(ps = ps1_rela, Micromet = Micromet, dist = dist)

  } else {
    pairResult = "no result"
  }


  #----绘图#----
  map = as.data.frame(phyloseq::sample_data(ps1_rela))
  map$Group = as.factor(map$Group)
  colbar = length(levels(map$Group))

  points = cbind(points, map[match(rownames(points), rownames(map)), ])
  # head(points)
  points$ID = row.names(points)


  if (method %in% c("DCA", "CCA", "RDA",  "MDS", "PCoA","PCA","LDA")) {
    p2 =ggplot(points, aes(x=x, y=y, fill = Group)) +
      geom_point(alpha=.7, size=5, pch = 21) +
      labs(x=paste0(method," 1 (",format(100*eig[1]/sum(eig),digits=4),"%)"),
           y=paste0(method," 2 (",format(100*eig[2]/sum(eig),digits=4),"%)"),
           title=title1) +
      stat_ellipse(linetype=2,level=0.68,aes(group=Group, colour=Group))

    p3 = p2+ggrepel::geom_text_repel(aes(label=points$ID),size = 5)
    p3
  }


  if (method %in% c("NMDS")) {
    p2 =ggplot(points, aes(x=x, y=y, fill = Group)) +
      geom_point(alpha=.7, size=5, pch = 21) +
      labs(x=paste(method,"1", sep=""),
           y=paste(method,"2",sep=""),
           title=stress)+
      stat_ellipse( linetype = 2,level = 0.65,aes(group  =Group, colour =  Group))

    p3 = p2+ggrepel::geom_text_repel( aes(label=points$ID),size=4)
    p3
    if (method %in% c("t-sne")) {
      supp_lab = labs(x=paste(method,"1", sep=""),y=paste(method,"2",sep=""),title=title)
      p2 = p2 + supp_lab
      p3 = p3 + supp_lab
    }
    eig = NULL
  }

  if (method %in% c("t-sne")) {
    p2 =ggplot(points, aes(x=x, y=y, fill = Group)) +
      geom_point(alpha=.7, size=5, pch = 21) +
      labs(x=paste(method,"1", sep=""),
           y=paste(method,"2",sep=""))+
      stat_ellipse( linetype = 2,level = 0.65,aes(group  =Group, colour =  Group))

    p3 = p2+ggrepel::geom_text_repel( aes(label=points$ID),size=4)
    p3
    if (method %in% c("t-sne")) {
      supp_lab = labs(x=paste(method,"1", sep=""),y=paste(method,"2",sep=""),title=title1)
      p2 = p2 + supp_lab
      p3 = p3 + supp_lab
    }
    eig = NULL
  }

  # 返回结果：标准图，数据，标签图，成对比较结果，整体结果
  return(list(p2,points,p3,pairResult,title1,eig))
}
