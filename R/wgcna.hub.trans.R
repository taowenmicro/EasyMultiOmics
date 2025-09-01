#' @title Perform WGCNA analysis on transcriptome functional composition table
#' @description This function conducts various steps of WGCNA analysis including data preprocessing, module identification,
#' trait correlation analysis, network visualization, and gene selection based on module-trait relationships.
#' @param ps A phyloseq format file used as an alternative for the input containing transcriptome
#' functional composition table, tax, and sample metadata.
#' @param data_type Transcriptome data type,default is "TPM".
#' If the input data is counts,log2(counts+1) will be executed;
#' If the input data is not counts, e.g., TPM、fpkm, there will be not convert again
#' @param group A string specifying the grouping variable in the sample metadata. Default is `"Group"`.
#' @param WGCNApath Path to save the results and plots from the analysis process,.
#' Default is `"NULL"`,this means that the the results and plots of the analysis process don't need to xported.
#' If you need to export, please modify it to the specified file output path.
#' @param alltraits Traits data, a data frame or matrix,rows are different samples,cols are different traits.
#' Default is `"NULL"`,this means that the grouping variable in the sample metadata will be considered the trait that needs attention.
#' If you have traits data , please modify it to the specified name of your traits data.
#' @param sample_filter Logical, whether to remove abnormal outlier samples,Default is `"TRUE"`
#' @param clustsample_thresold Threshold for sample clustering. This parameter only takes effect when sample_filter is true
#' @param minModuleSize Minimum module size for module identification.Default is `"30"`.
#' @param MEDissThres Similarity threshold for merging modules.Modules with similarity exceeding this threshold will be merged into one module.Default is `"0.75"`.
#' @param topnSelect_module  Number of top modules to select, based on the absolute correlation between traits and modules.Default is `"5"`.
#' @param nSelect_genes Number of genes to select for heatmap.Default is `"400"`.
#' @return A list containing the top module-trait relationships and selected genes information based on MM and GS.
#' \item {moduleTraitCor_top}{Data frame containing the selected top modules with all traits.}
#' \item {xx3}{Data frame containing important genes by MM and GS in the selected top modules with all traits.}
#' @export
#' @examples
#' \dontrun{
#' library(WGCNA)
#' library(dplyr)
#' res <- wgcna.trans(ps=ps.trans, data_type="TPM",group="Group",WGCNApath=NULL,
#'                    alltraits=NULL,sample_filter=TRUE,clustsample_thresold=50000,
#'                    minModuleSize = 30,MEDissThres = 0.25,topnSelect_module=5,
#'                    nSelect_genes = 400)
#' head(res[[1]])
#' head(res[[2]])
#' }
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njau.edu.cn}
wgcna.hub.trans <- function(ps=ps.trans, data_type="TPM",group="Group",WGCNApath=NULL,
                            alltraits=NULL,sample_filter=FALSE,clustsample_thresold=50000,
                            minModuleSize = 30,MEDissThres = 0.25,topnSelect_module=5,
                            nSelect_genes = 400){
  if (!is.null(WGCNApath)) {a = Sys.Date() %>% as.character()
  a = gsub("-",".",a)
  WGCNApath = paste(WGCNApath,"/",a,sep = "");WGCNApath
  dir.create(WGCNApath,recursive = T)}
  options(stringsAsFactors = FALSE)  #开启多线程
  femData = data.frame(phyloseq::otu_table(ps.trans)) #载入基因表达量数据
  if (is.null(alltraits)) {
    groupdata=data.frame(phyloseq::sample_data(ps.trans))
    alltraits <- data.frame(row.names = rownames(groupdata))
    for (i in unique(groupdata[,group])) {
      alltraits <- alltraits%>%dplyr::mutate( !!sym(i):=dplyr::case_when(groupdata[,group] == i ~ 1, TRUE ~ 0))
    }
  }
  #要转化为行为样本，列为基因，
  if (data_type=="counts") {
    datExpr0 <- t(log2(femData+1))
  }else{
    datExpr0 <- t(femData)
  }
  #1.2 检查基因缺失值和识别离群值（异常值）----
  gsg = goodSamplesGenes(datExpr0, verbose = 3);
  gsg$allOK
  if (!gsg$allOK){
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(colnames(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
  }
  #1.3聚类所有样本，观察是否有样本离群值或异常值----
  sampleTree = hclust(dist(datExpr0), method = "average")
  if (!is.null(WGCNApath)) {
    pdf(file = paste(WGCNApath,"/","1_sampleClustering.pdf",sep = ""), width = 12, height = 9)
  }
  sizeGrWindow(12,9) #视图
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  if (sample_filter==TRUE){abline(h = clustsample_thresold, col = "red")}
  dev.off()
  #裁剪离群样本
  if (sample_filter==TRUE){
    clust = cutreeStatic(sampleTree, cutHeight = clustsample_thresold, minSize = 10)
    table(clust)
    keepSamples = (clust==1)
    datExpr0 = datExpr0[keepSamples, ]
    sampleTree = hclust(dist(datExpr0), method = "average")
    plot(sampleTree)}
  #2不同样本性状变化----
  traits <- alltraits[rownames(datExpr0),]
  traitColors = numbers2colors(traits, signed = FALSE)
  if (!is.null(WGCNApath)) {
    pdf(file=paste(WGCNApath,"/","2_Sample dendrogram and trait heatmap.pdf",sep = ""),width=12,height=12)}
  plotDendroAndColors(sampleTree, traitColors,
                      groupLabels = names(traits),
                      main = "Sample dendrogram and trait heatmap")
  dev.off()
  #3构建表达网络----
  enableWGCNAThreads(nThreads = 6)  #开启多线程
  # Choose a set of soft-thresholding powers幂指数计算范围，
  powers = c(1:30)
  # Call the network topology analysis function估计最佳power值，自动估计和计算
  sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
  if (!is.null(WGCNApath)) {
    pdf(file=paste(WGCNApath,"/","3_Scale independence.pdf",sep = ""),width=9,height=5)}
  par(mfrow = c(1,2))
  cex1 = 0.9
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
  ######chose the softPower，估计出最佳的power值，转化为邻接矩阵
  softPower =sft$powerEstimate
  # softPower=7
  softPower
  #转化为邻接矩阵,
  adjacency = adjacency(datExpr0, power = softPower)

  ##### Turn adjacency into topological overlap
  TOM = TOMsimilarity(adjacency);
  #TOM相异矩阵
  dissTOM = 1-TOM
  # Call the hierarchical clustering function将TOM相异性矩阵进行聚类，相异性表示距离看相关
  geneTree = hclust(as.dist(dissTOM), method = "average");
  # Plot the resulting clustering tree (dendrogram)

  #sizeGrWindow(12,9)
  if (!is.null(WGCNApath)) {
    pdf(file=paste(WGCNApath,"/","4_Micro clustering on TOM-based dissimilarity.pdf",sep = ""),width=15,height=11)}
  plot(geneTree, xlab="", sub="", main = "Micro clustering on TOM-based dissimilarity",
       labels = FALSE, hang = 0.04)
  dev.off()

  # Module identification using dynamic tree cut:计算模块
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize);
  #，0代表什么模块都不是颜色就是灰色，显示哪个模块有多少个基因
  table(dynamicMods)
  #去除灰色
  # Convert numeric lables into colors，匹配颜色
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  # Plot the dendrogram and colors underneath
  if (!is.null(WGCNApath)) {
    pdf(file=paste(WGCNApath,"/","5_Dynamic Tree Cut.pdf",sep = ""),width=12,height=9)}
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  dev.off()

  #小的模块相似性比较高，需要进行合并，这里引入模块特征基因和模块特征向量
  # Calculate eigengenes模块合并，首先计算模块特征向量
  MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
  MEs = MEList$eigengenes
  MEs <- MEs[,which(!colnames(MEs)=="MEgrey")]
  which(is.na(MEs)==T)
  # Calculate dissimilarity of module eigengenes，转化相异性
  MEDiss = 1-cor(MEs);
  # Cluster module eigengenes，聚类并绘制聚类树
  METree = hclust(as.dist(MEDiss), method = "average")
  # Plot the result
  #sizeGrWindow(7, 6)
  if (!is.null(WGCNApath)) {
    pdf(file=paste(WGCNApath,"/","6_Clustering of module eigengenes.pdf",sep = ""),width=7,height=6)}
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  # Plot the cut line into the dendrogram，
  abline(h=MEDissThres, col = "red")
  dev.off()


  # Call an automatic merging function合并模块命令cutHeight是设置合并高度
  #对应相似性即为0.75，按照相似性为0.75的模块合并为一个模块
  merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors
  unique(mergedColors)
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs

  #sizeGrWindow(12, 9)
  #合并之后模块的展示，每个颜色合并成一个模块，合并前和合并后一起展示
  if (!is.null(WGCNApath)) {
    pdf(file=paste(WGCNApath,"/","7_merged dynamic.pdf",sep = ""), width = 9, height = 6)}
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()

  # Rename to moduleColors，查看一下合并后有多少个模块，每个模块有多少基因，灰色grey的代表这些基因不属于任何一个模块
  moduleColors = mergedColors
  table(moduleColors)
  # Construct numerical labels corresponding to the colors将模块排序，将编号进行对应,将灰色定义为0
  colorOrder = c("grey", standardColors(100))
  moduleLabels = match(moduleColors, colorOrder)-1
  table(moduleLabels)
  MEs = mergedMEs
  #这部分只要保存就可以了，以后不需要再次运行以上的代码
  # Save module colors and labels for use in subsequent parts，网络相关数据进行存储
  if (!is.null(WGCNApath)) {
    save(MEs, TOM, dissTOM,  moduleLabels, moduleColors, geneTree, sft, file = paste(WGCNApath,"/","networkConstruction-stepByStep.RData",sep = ""))}

  # load(paste(WGCNApath,"/","networkConstruction-stepByStep.RData",sep = ""))

  #8比较模块（特征模块向量）和性状的相关分析#--------
  MEs <- MEs[,which(!colnames(MEs)=="MEgrey")]
  nGenes = ncol(datExpr0)
  nSamples = nrow(datExpr0)
  moduleTraitCor = cor(MEs, traits, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  #可视化
  #sizeGrWindow(10,6)
  if (!is.null(WGCNApath)) {
    pdf(file=paste(WGCNApath,"/","8_Module-trait relationships.pdf",sep = ""),width=7,height=15)}
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "")

  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3))
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(traits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(traits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  ##8.3输出选定重要模块----
  moduleTraitCor_top <- data.frame(matrix(ncol =topnSelect_module, nrow =ncol(moduleTraitCor)))
  rownames(moduleTraitCor_top) <- colnames(moduleTraitCor)
  for (i in colnames(moduleTraitCor)) {
    moduleTraitCor_top[i,] <- names(moduleTraitCor[,i][
      order(abs(moduleTraitCor[,i]),decreasing = TRUE)[1:topnSelect_module]])
  }
  colnames(moduleTraitCor_top) <- paste("top",1:topnSelect_module,sep = "_")
  #9MM和GS确定重要基因----
  ##9.1-模块和基因相关性----
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) = paste("MM", modNames, sep="")
  names(MMPvalue) = paste("p.MM", modNames, sep="")
  ##9.2-性状和基因相关性----
  traitNames=names(traits)
  geneTraitSignificance = as.data.frame(cor(datExpr0, traits, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
  names(GSPvalue) = paste("p.GS.", traitNames, sep="")
  ##9.3-MM和GS----
  if (!is.null(WGCNApath)) {
    WGCNApath1 <- paste(WGCNApath,"/9_性状与模块",sep="")
    dir.create(WGCNApath1)}
  # trait <- "NC"
  # module <- "skyblue"
  for (trait in traitNames){
    traitColumn=match(trait,traitNames)
    for (module in modNames){
      column = match(module, modNames)
      moduleGenes = moduleColors==module
      xx <- data.frame(row.names =rownames(geneModuleMembership)[moduleGenes] ,
                       MM=geneModuleMembership[moduleGenes, column],
                       GS=geneTraitSignificance[moduleGenes, traitColumn])
      xx1 <- xx[which(abs(xx$MM)>0.8&abs(xx$GS)>0.5),]
      xx1 <- xx1[order(xx1$GS,decreasing = T),]
      if (!is.null(WGCNApath)) {
        hubpath <- paste(WGCNApath1,"/",trait, "_", module,"_hub.csv",sep="")
        write.csv(xx1, hubpath,quote = F)}
      if (nrow(geneModuleMembership[moduleGenes,]) > 1){
        if (!is.null(WGCNApath)) {
          outPdf=paste(WGCNApath1,"/",trait, "_", module,".pdf",sep="")
          pdf(file=outPdf,width=7,height=7)}
        par(mfrow = c(1,1))
        verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                           abs(geneTraitSignificance[moduleGenes, traitColumn]),
                           xlab = paste("Module Membership in", module, "module"),
                           ylab = paste("Gene significance for ",trait),
                           main = paste("Module membership vs. gene significance\n"),
                           cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
        abline(v=0.8,h=0.5,col="red")
        dev.off()
      }
    }
  }
  ##9.4导出top模块重要基因----
  for (trait in rownames(moduleTraitCor_top)){
    traitColumn=match(trait,traitNames)
    for (module in substring(unique(moduleTraitCor_top[trait,]),3)){
      column = match(module, modNames)
      moduleGenes = moduleColors==module
      xx <- data.frame(row.names =rownames(geneModuleMembership)[moduleGenes] ,
                       MM=geneModuleMembership[moduleGenes, column],
                       GS=geneTraitSignificance[moduleGenes, traitColumn])
      xx1 <- xx[which(abs(xx$MM)>0.8&abs(xx$GS)>0.5),]
      if (nrow(xx1) == 0) {
        next}
      xx1 <- xx1[order(xx1$GS,decreasing = T),]
      xx1$type <- paste(trait,"top",which(substring(moduleTraitCor_top[trait,],3)==module),
                        module,sep = "_")
      if (module==substring(unique(moduleTraitCor_top[trait,]),3)[1]) {
        xx2 <- xx1}else{
          xx2 <- rbind(xx2,xx1)}}
    if (nrow(xx2) == 0) {
      next}
    if (trait==rownames(moduleTraitCor_top)[1]) {
      xx3 <- xx2}else{
        xx3 <- rbind(xx3,xx2)}}
  #10-导出每个模块的网络图文件----
  if (!is.null(WGCNApath)) {
    cytoDir=paste(WGCNApath,"/10_Cytoscape/",sep="")
    dir.create(cytoDir)
    for (mod in 1:nrow(table(moduleColors))){ #【对每个模块进行循环】
      modules = names(table(moduleColors))[mod]
      probes = colnames(datExpr0)
      inModule = (moduleColors == modules)
      modProbes = probes[inModule]
      modGenes = modProbes
      modTOM = TOM[inModule, inModule]
      dimnames(modTOM) = list(modProbes, modProbes)
      edges_File = paste("CytoscapeInput-edges-", modules , ".txt", sep="")
      nodes_File = paste("CytoscapeInput-nodes-", modules, ".txt", sep="")
      outEdge=paste(cytoDir,edges_File,sep="\\")
      outNode=paste(cytoDir,nodes_File,sep="\\")
      cyt = exportNetworkToCytoscape(modTOM,
                                     edgeFile = outEdge,
                                     nodeFile = outNode,
                                     weighted = TRUE,
                                     threshold = 0.02,
                                     nodeNames = modProbes,
                                     altNodeNames = modGenes,
                                     nodeAttr = moduleColors[inModule])
    }}
  #11-导出每个模块所包含的基因----
  if (!is.null(WGCNApath)) {
    moduleDir=paste(WGCNApath,"/11_模块基因/",sep="")
    dir.create(moduleDir)}
  for (mod in 1:nrow(table(moduleColors))){
    modules = names(table(moduleColors))[mod]
    probes = colnames(datExpr0)
    inModule = (moduleColors == modules)
    modGenes = probes[inModule]
    if (!is.null(WGCNApath)) {
      write.table(modGenes, file =paste0(moduleDir,modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
    }}
  #--12-基因相关热图----
  # Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
  plotTOM = dissTOM^7
  # Set diagonal to NA for a nicer plot
  diag(plotTOM) = NA

  if (nSelect_genes=="All") {pdf(file= paste(WGCNApath,"/","12_Network heatmap plot_all gene.pdf",sep = ""),width=9, height=9)
    TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
    dev.off()
  }else{set.seed(10)
    select = sample(nGenes, size = nSelect_genes)
    selectTOM = dissTOM[select, select]
    # There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
    selectTree = hclust(as.dist(selectTOM), method = "average")
    selectColors = moduleColors[select]
    # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
    # the color palette; setting the diagonal to NA also improves the clarity of the plot
    plotDiss = selectTOM^7
    diag(plotDiss) = NA
    if (!is.null(WGCNApath)) {
      pdf(file=paste(WGCNApath,"/","12_Network heatmap plot_selected genes.pdf",sep = ""),width=9, height=9)}
    TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
    dev.off()}
  return(list(moduleTraitCor_top,xx3))
}
