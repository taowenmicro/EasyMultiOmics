#' @title Generate a Circular Microbial Network
#'
#' @description
#' This function generates a circular microbial network based on microbial abundance data in a `phyloseq` object.
#' It uses the `ggClusterNet` package to compute relationships between operational taxonomic units (OTUs)
#' and visualizes the network with sample-level metadata.
#' @param ps A `phyloseq` object containing microbiome data.
#' @param N Numeric. A threshold for selecting nodes in the network. Default is `0.5`.
#' @param fill A character string specifying the taxonomic rank used to color the nodes. Default is `"Phylum"`.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{plot}{A `ggplot2` object representing the microbial network.}
#'   \item{plotdata}{A data frame containing node coordinates and associated metadata.}
#' }
#'
#' @details
#' This function builds a microbial network using a circular layout. OTUs are placed hierarchically based on their
#' abundance and relationships with samples. The edges represent connections between OTUs and sample groups,
#' while nodes represent microbial taxa or samples.
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' The function includes several steps:
#' \itemize{
#'   \item Calculates the microbial relationships using `ggClusterNet::div_network`.
#'   \item Uses a progressive circle packing algorithm to position OTUs and groups.
#'   \item Combines abundance data and taxonomic annotations for visualization.
#'   \item Generates a `ggplot2` plot showing the network with colored nodes and edges.
#' }
#' @examples
#' result = ven.network.micro(ps = ps.16s,N = 0.5,fill = "Phylum")
#' p14  = result[[1]]
#' p14
#' dat = result[[2]]
#' @export

ven.network.micro = function(
  ps =ps,
  N = 0.5,
  fill = "Phylum"
  ){
  result = ggClusterNet::div_network(ps)

  edge = result[[1]]
  data = result[[3]]
  print("1")
  # result <- ggClusterNet::div_culculate(table = result[[3]],
  #                                       distance = 1.1,distance2 = 1.5,distance3 = 1.3,order = FALSE)
  # result <- div_culculate(table = result[[3]],distance = 1,distance2 = 1.2,distance3 = 1.1,order = FALSE)
  table = result[[3]]
  distance = 1.1
  distance2 = 1.5
  distance3 = 1.3
  order = FALSE


  all = table[rowSums(table) == dim(table)[2],]
  allnum = dim(table[rowSums(table) == dim(table)[2],])[1]
  N = allnum
  if (N == 0) {
    print("all cover was 0,so,can not continue")
  }else {
    if (order== TRUE) {
      packing <- packcircles::circleProgressiveLayout(rep(1,N))

    }
    if (order== FALSE) {
      packing <- packcircles::circleProgressiveLayout(runif( min = 1, max = 10,n=N))

    }

    data <- packcircles::circleLayoutVertices(packing)  %>% dplyr::group_by(id) %>%
      dplyr::summarise(x = mean(x),y = mean(y))  %>%
      dplyr::select(-id)  %>%
      as.data.frame()
    dim(data)
    r0 = max(data$x) - min(data$x)
    row.names(data ) = row.names(all)
    data$elements = row.names(data )
    colnames(data)[1:2] = c("X1","X2")
    allxy = data


    #---culculate da
    r = distance*r0
    #--Calculate angle according to group
    arg = seq(0,360,360/(dim(table)[2]))
    x= rep(0,dim(table)[2])
    y = rep(0,dim(table)[2])
    for (i in 1:dim(table)[2]) {

      x[i] = r* sin(arg[i]* 3.14/180)

      y[i] = r* cos(arg[i]* 3.14/180)
    }

    da0 = data.frame(x = x,y = y)

    #---------
    r = distance2*r0
    x= rep(0,dim(table)[2])
    y = rep(0,dim(table)[2])
    for (i in 1:dim(table)[2]) {
      x[i] = r* sin(arg[i]* 3.14/180)
      y[i] = r* cos(arg[i]* 3.14/180)
    }
    da1 = data.frame(x = x,y = y)
    da1
    # i = 1
    for (i in 1:dim(table)[2]) {

      table = table[rowSums(table) != dim(table)[2],]
      N = length(table[,i][table[,i] != 0])

      if (N != 0) {
        if (order== TRUE) {
          packing <- packcircles::circleProgressiveLayout(rep(1,N))
        } else {
          packing <- packcircles::circleProgressiveLayout(runif(min = 1, max = 10,n=N))
        }
        data <- packcircles::circleLayoutVertices(packing)  %>% dplyr::group_by(id) %>%
          dplyr::summarise(x = mean(x),y = mean(y))  %>%
          dplyr::select(-id)  %>%
          as.data.frame()
        data <- data.frame(x = data$x + da1[i,1] * distance3,y = data$y + da1[i,2]* distance3)
        row.names(data ) = row.names(table[table[,i] != 0,])
        data$elements = row.names(data )
        colnames(data)[1:2] = c("X1","X2")
      } else{
        data = NULL
      }

      if (i == 1) {
        oridata = data
      }
      if (i != 1) {
        oridata = rbind(oridata,data)
      }
    }

    oridata = rbind(oridata,allxy)
    head(oridata )

    rownames(da0) = colnames(table)
    da0$elements = row.names(da0 )
    colnames(da0)[1:2] = c("X1","X2")
    # plotdata = rbind(oridata,da0)
    plotdata=  oridata

    edg <-plotdata %>%
      dplyr::rename(X1 = X1,Y1 = X2,OTU = elements) %>%
      # select(X) %>%
      dplyr::right_join(edge,by = c("OTU"="source"))

    edge2 <-da0 %>%
      dplyr::rename(X2 = X1,Y2 = X2,Group = elements) %>%
      dplyr::right_join(edg,by = c("Group"="target"))

    result  =  list(edge2,plotdata,da0)
  }

   edge = result[[1]]

  plotdata = result[[2]]

  #--这部分数据是样本点数据
  groupdata <- result[[3]]
  # table(plotdata$elements)
  node =  plotdata[plotdata$elements == unique(plotdata$elements), ]

  otu_table = as.data.frame(t(ggClusterNet::vegan_otu(ps)))
  tax_table = as.data.frame(ggClusterNet::vegan_tax(ps))
  res = merge(node,tax_table,by = "row.names",all = F)
  row.names(res) = res$Row.names
  res$Row.names = NULL
  plotcord = res
  xx = data.frame(mean  =rowMeans(otu_table))
  plotcord = merge(plotcord,xx,by = "row.names",all = FALSE)
  head(plotcord)
  row.names(plotcord) = plotcord$Row.names
  plotcord$Row.names = NULL
  p = ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2),
                              data = edge, size = 0.3,color = "yellow") +
    geom_point(aes(X1, X2,fill = !!sym(fill),size =mean ),pch = 21, data = plotcord) +
    geom_point(aes(X1, X2),pch = 21, data = groupdata,size = 5,fill = "blue",color = "black") +
    geom_text(aes(X1, X2,label = elements ), data = groupdata,hjust = 1,vjust = -1) +
    theme_void()
  p


  return(list(plot = p,plotdata = plotcord))
}













