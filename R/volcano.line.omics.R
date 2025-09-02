
#' @title Volcano Plot with Correlation and Random Forest Analysis for Omics Data
#'
#' @description
#' This function compares two omics datasets (`ps1` and `ps2`) by generating volcano plots that depict the differential
#' expression based on log fold changes and significance (p-values). It incorporates statistical analysis (e.g., Pearson
#' correlation) and visualizes the relationships between the datasets. The plots color points to represent enriched,
#' depleted, or non-significant variables, and uses edge networks based on correlation values between variables.
#' Additionally, the function performs Random Forest analysis to evaluate the importance of environmental or experimental
#' variables influencing the data.
#'
#' @examples
#' \dontrun{
#' # Assuming 'ps.trans' and 'ps.16s' are valid phyloseq objects with data
#' ps1 = ps.trans %>% filter_OTU_ps(500)
#' ps2 = ps.16s %>% filter_OTU_ps(500)
#'
#' res = volcano.line.omics(ps1 = ps1, ps2 = ps2,
#'                           lab.1 = "tran", lab.2 = "16s",
#'                           group = "Group", r.threshold = 0.8,
#'                           p.threshold = 0.05, col.point = NULL,
#'                           col.line = NULL, method.cor = "pearson",
#'                           n = 50)
#'
#' # View the volcano plot results
#' res[[1]]$WT_OE
#' res[[1]]$WT_KO
#' res[[1]]$OE_KO
#' }
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}

volcano.line.omics = function(
    ps1 = ps1,
    ps2 = ps2,
    lab.1 = "tran",
    lab.2 = "16s",
    group = "Group",
    r.threshold = 0.8,
    p.threshold= 0.05,
    col.point = NULL,
    col.line = NULL,
    method.cor = "pearson",
    n = 50
){


  res1 = EdgerSuper.trans(ps = ps1,group  =  group,artGroup = NULL)
  head(res1)
  map = sampleData(ps.trans)
  ids = map$Group %>%
    unique() %>%
    combn(2)
  #  两两组合的组做差异，展示火山图也是展示两两组合

  ids[,i]
  plot1 = list()
  plot2 = list()
  plot3 = list()
  dat1 = list()
  dat2 = list()
  dat3 = list()
  for (i in 1:dim(ids)[2]) {
    dat = res1[[2]]
    head(dat)
    #  计计算OTU之间相关
    n1 = cor_Big_micro(ps1 %>% subset_samples.wt(group, ids[,i]) %>% remove.zero(),
                       r.threshold = r.threshold,
                       method = method.cor
    )

    cor = n1[[1]]
    a = colnames(dat)[str_detect(colnames(dat),"logFC")]
    b = colnames(dat)[str_detect(colnames(dat),ids[2,i])]
    c = colnames(dat)[str_detect(colnames(dat),ids[1,i])]
    d = intersect(a,b)
    e = intersect(d,c)

    a = colnames(dat)[str_detect(colnames(dat),"p")]
    b = colnames(dat)[str_detect(colnames(dat),ids[2,i])]
    c = colnames(dat)[str_detect(colnames(dat),ids[1,i])]
    d = intersect(a,b)
    e2 = intersect(d,c)

    a = colnames(dat)[str_detect(colnames(dat),"level")]
    b = colnames(dat)[str_detect(colnames(dat),ids[2,i])]
    c = colnames(dat)[str_detect(colnames(dat),ids[1,i])]
    d = intersect(a,b)
    e3 = intersect(d,c)

    node = data.frame(elements = row.names(dat),
                      X1 = dat[,e],
                      X2 =dat[,e2])
    node$X2 = -log2(node$X2)
    node0 = node
    node$X1[node$X1>0] = - node$X1[node$X1>0]


    node2 = data.frame(elements =row.names(dat),
                       X1 = dat[,e],
                       X2 =dat[,e2],
                       class = dat[,e3]
    )
    node20 = node2
    node2$X2 = -log2(node2$X2)
    node20 = node2
    node2$X1[node2$X1>0] = - node2$X1[node2$X1>0]
    node2.1 = node2

    edge = edgeBuild(cor = cor,node = node0)
    head(edge)



    if (is.null(col.point)&is.null(col.line)) {
      col1 = c("red","blue")
      names(col1) = c("+","-")
      cols = c('blue2','red2', 'gray30')
      names(cols) = c("depleted","enriched","nosig")
    } else{
      col1 = col.point
      cols = col.line
    }

    p1 <- ggplot() +
      geom_point(data = node20,aes(x =X1 ,y = X2, colour=class),size = 3,pch = 21) +
      scale_color_manual(values = cols) +
      ggnewscale::new_scale_color() +

      geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                   data = edge, size = 0.5) +
      geom_hline(yintercept=2.5,
                 linetype=4,
                 color = 'black',
                 size = 0.5) +
      geom_vline(xintercept=c(-1,1),
                 linetype=3,
                 color = 'black',
                 size = 0.5) +
      # ggrepel::geom_text_repel( aes(x =logFC ,y = -log2(p), label=Genus), size=1) +
      scale_color_manual(values = col1)+
      theme_bw()
    # p1
    plot1[[i]] = p1
    names(plot1)[i] = paste(ids[1,i],ids[2,i],sep = "_")
    dat1[[i]] = list(node = node20,edge = edge)
    names(dat1)[i] = paste(ids[1,i],ids[2,i],sep = "_")


    res = EdgerSuper.trans(ps = ps2,group  = group,artGroup = NULL)
    head(res)
    dat = res[[2]]
    head(dat)

    # map = sampleData(ps2)
    # ids = map$Group %>%
    #   unique() %>%
    #   combn(2)
    # i = 1
    # ids[,i]


    n1 = cor_Big_micro(ps2 %>% subset_samples.wt(group, ids[,i]) %>% remove.zero(),
                       r.threshold = 0.85
    )
    cor = n1[[1]]

    a = colnames(dat)[str_detect(colnames(dat),"logFC")]
    b = colnames(dat)[str_detect(colnames(dat),ids[2,i])]
    c = colnames(dat)[str_detect(colnames(dat),ids[1,i])]
    d = intersect(a,b)
    e = intersect(d,c)

    a = colnames(dat)[str_detect(colnames(dat),"p")]
    b = colnames(dat)[str_detect(colnames(dat),ids[2,i])]
    c = colnames(dat)[str_detect(colnames(dat),ids[1,i])]
    d = intersect(a,b)
    e2 = intersect(d,c)

    a = colnames(dat)[str_detect(colnames(dat),"level")]
    b = colnames(dat)[str_detect(colnames(dat),ids[2,i])]
    c = colnames(dat)[str_detect(colnames(dat),ids[1,i])]
    d = intersect(a,b)
    e3 = intersect(d,c)

    node = data.frame(elements =row.names(dat),
                      X1 = dat[,e],
                      X2 =dat[,e2])
    node$X2 = -log2(node$X2)
    head(node)
    node0 = node
    node$X1[node$X1<0] = - node$X1[node$X1<0]


    node2 = data.frame(elements =row.names(dat),
                       X1 = dat[,e],
                       X2 =dat[,e2],
                       class = dat[,e3]
    )
    node2$X2 = -log2(node2$X2)
    node20 = node2

    node2$X1[node2$X1<0] = - node2$X1[node2$X1<0]
    node2.2 = node2
    edge = edgeBuild(cor = cor,node = node0)
    head(edge)


    # col1 = c("red","blue")
    # names(col1) = c("+","-")
    # cols = c('blue2','red2', 'gray30')
    # names(cols) = c("depleted","enriched","nosig")
    head(node2)

    p2 <- ggplot() +
      geom_point(data = node20,aes(x =X1 ,y = X2, colour=class),size = 3,pch = 21) +
      scale_color_manual(values = cols) +
      ggnewscale::new_scale_color() +

      geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                   data = edge, size = 0.5) +
      geom_hline(yintercept=2.5,
                 linetype=4,
                 color = 'black',
                 size = 0.5) +
      geom_vline(xintercept=c(-1,1),
                 linetype=3,
                 color = 'black',
                 size = 0.5) +
      # ggrepel::geom_text_repel( aes(x =logFC ,y = -log2(p), label=Genus), size=1) +
      scale_color_manual(values = col1)+
      theme_bw()

    plot2[[i]] = p2
    names(plot2)[i] = paste(ids[1,i],ids[2,i],sep = "_")
    dat2[[i]] = list(node = node20,edge = edge)
    names(dat2)[i] = paste(ids[1,i],ids[2,i],sep = "_")



    if (is.null(lab.1)&is.null(lab.2)) {
      node2.1$gro = "omics1"
      node2.2$gro = "omics1"
    } else{
      node2.1$gro = lab.1
      node2.2$gro = lab.2
    }

    node3 = rbind(node2.1,node2.2)
    head(node3)

    if (length(node3$X1[node3$X1 >5])> 0) {
      node3$X1[node3$X1 >5] = node3$X1[node3$X1 >5] %>% log2()
    }

    if (length(node3$X1[node3$X1 < -5])> 0) {
      node3$X1[node3$X1 < -5] = -(-node3$X1[node3$X1 < -5] %>% log2())
    }



    head(node3)
    node3$elements2 = paste(node3$gro,node3$elements,sep = "_")
    node3$elements = node3$elements2



    ps.all = merge.ps(ps1 = ps1,
                      ps2 = ps2,
                      N1 = 0,
                      N2 = 0,
                      scale = TRUE,
                      onlygroup = TRUE,#不进行列合并，只用于区分不同域
                      dat1.lab = lab.1,
                      dat2.lab = lab.2)


    #  计计算OTU之间相关
    result <- corBiostripeBig(ps = ps.all,
                              r.threshold=r.threshold,
                              p.threshold=p.threshold,
                              method = method.cor)

    # extract cor matrix
    cor = result[[1]]
    # n1 = cor_Big_micro(ps.all,
    #                    r.threshold = 0.95
    # )
    #
    # cor = n1[[1]]

    edge = edgeBuild(cor = cor,node = node3)
    head(edge)
    dim(edge)

    # library(ggnewscale)
    col1 = c("#9E0142", "#3288BD")
    names(col1) = c("+","-")
    cols = c("#66C2A5","#F46D43", 'gray30')
    names(cols) = c("depleted","enriched","nosig")
    head(node3)
    head(edge)
    colnames(edge)[length(colnames(edge))] = "cor.col"


    labs = node3 %>% dplyr::filter(class != "nosig") %>% arrange(desc(X2)) %>% head(n)


    p3 <- ggplot() +
      geom_jitter(data = node3,aes(x =X1 ,y = X2, colour=class,shape = gro),size = 3) +
      scale_color_manual(values = cols) +
      ggrepel::geom_text_repel(aes(x =X1 ,y = X2,label = elements2),data = labs)+
      ggnewscale::new_scale_color() +

      geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor.col)),
                   data = edge, size = 0.5,alpha = 0.3) +
      geom_hline(yintercept=2.5,
                 linetype=4,
                 color = 'black',
                 size = 0.5) +
      geom_vline(xintercept=c(-1,1),
                 linetype=3,
                 color = 'black',
                 size = 0.5) +
      # ggrepel::geom_text_repel( aes(x =logFC ,y = -log2(p), label=Genus), size=1) +
      scale_color_manual(values = col1)+
      theme_bw()
    p3
    plot3[[i]] = p3
    names(plot3)[i] = paste(ids[1,i],ids[2,i],sep = "_")
    dat3[[i]] = list(node = node20,edge = edge)
    names(dat3)[i] = paste(ids[1,i],ids[2,i],sep = "_")

  }


  return(list(plot3,plots = list(plot1,plot2,plotdata = list(dat1,dat2,dat3))))

}

