

# #--维恩网络#-------
# library(ggClusterNet)
# library(phyloseq)
# library(ggrepel)
#
# biospath = paste(otupath,"/biospr_network_Ven/",sep = "")
# dir.create(biospath)
#
# result = ven.network(
#     ps = ps,
#     N = 0.5,
#     fill = "Phylum"
#     )
#
# p  = result[[1]]
#
# data = result[[2]]
#
# filename = paste(biospath,"/","biostr_Ven_network.species.several.pdf",sep = "")
# ggsave(filename,p,width = (15),height = (12))
# filename = paste(biospath,"/","biostr_Ven_network.jpg",sep = "")
# ggsave(filename,p,width = (15),height = (12))
#
# filename = paste(biospath,"Ven.network.all.csv",sep = "")
# write.csv(data,filename)
#
#
# detach("package:ggClusterNet")
# detach("package:phyloseq")


div_network2 = function (ps, group = "Group", flour = TRUE, N = 0.5) {
  mapping = as.data.frame(phyloseq::sample_data(ps))
  mapping = mapping[, group]
  colnames(mapping[, group]) <- "Group"
  sample_data(ps) = mapping
  ps_rela = phyloseq::transform_sample_counts(ps, function(x) x/sum(x))
  ps_rela
  aa = vegan_otu(ps)
  otu_table = as.data.frame(t(aa))
  count = aa
  countA = count
  dim(count)
  sub_design <- as.data.frame(phyloseq::sample_data(ps))
  count[count > 0] <- 1
  count2 = as.data.frame(count)
  iris.split <- split(count2, as.factor(sub_design$Group))
  iris.apply <- lapply(iris.split, function(x) colSums(x[]))
  iris.combine <- do.call(rbind, iris.apply)
  ven2 = t(iris.combine)
  for (i in 1:length(unique(phyloseq::sample_data(ps)$Group))) {
    aa <- as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,
                                                                1]
    bb = as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,
                                                               2]
    ven2[, aa] = ven2[, aa]/bb
  }
  ven2[ven2 < N] = 0
  ven2[ven2 >= N] = 1
  ven2 = as.data.frame(ven2)
  ven3 = as.list(ven2)
  ven2 = as.data.frame(ven2)
  if (flour == TRUE) {
    ven2 = ven2[rowSums(ven2) == dim(ven2)[2] | rowSums(ven2) ==
                  1, ]
  }
  tax_table = as.data.frame(vegan_tax(ps))
  otu_table = as.data.frame(t(vegan_otu(ps)))
  dim(otu_table)
  otu_net = merge(ven2, tax_table, by = "row.names", all = F)
  row.names(otu_net) = otu_net$Row.names
  otu_net$Row.names = NULL
  head(otu_net)
  OTU = as.matrix(otu_table)
  norm = t(t(OTU)/colSums(OTU))
  norm1 = norm %>% t() %>% as.data.frame()
  iris.split <- split(norm1, as.factor(mapping$Group))
  iris.apply <- lapply(iris.split, function(x) colSums(x))
  norm2 <- do.call(rbind, iris.apply) %>% t()
  head(norm2)
  colnames(norm2) = paste(colnames(norm2), "mean", sep = "")
  dim(otu_net)
  otu_net2 = merge(otu_net, norm2, by = "row.names", all = F)
  dim(otu_net2)
  head(otu_net2)
  colnames(otu_net2)[1] = "ID"
  sample_label = otu_net2[1:length(unique(mapping$Group)),
  ]
  sample_label$ID = unique(mapping$Group)
  point_label = rbind(otu_net2, sample_label)
  head(otu_net2)

  net_all = reshape2::melt(otu_net2, id = c("ID", rank_names(ps),
                                            paste(unique(mapping$Group), "mean", sep = "")))

  net_filter <- dplyr::filter(net_all, value != 0)

  net_fal = data.frame(source = net_filter$ID, target = net_filter$variable,
                       connect = rep("pp", nrow(net_filter)), value = net_filter$value,
                       Label = net_filter$ID)
  head(net_fal)
  return(list(net_fal, point_label, ven2))
}


ven.network.trans = function(
    ps =ps,
    N = 0.5,
    fill = "Level1"
){

  result = div_network2(ps)

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













