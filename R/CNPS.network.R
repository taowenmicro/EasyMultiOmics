#' @title Construct Co-occurrence Network for CNPS Data
#' @description This function constructs a co-occurrence network from CNPS gene
#' @param ps A phyloseq object containing OTU/ASV abundance data and taxonomy information
#' @param dat A data frame containing additional metadata or gene annotation information
#' @param id.0 Character string specifying the prefix for node IDs in the network (default: "C")
#'
#' @return A list containing network analysis results. Typically includes:
#' \itemize{
#'   \item Network graph object
#'   \item Node and edge attributes
#'   \item Network visualization
#' }
#' @examples
#' \dontrun{
#' # Example usage:
#' dat = db.cnps
#' network_result <- CNPS.network(
#'     ps = ps.trans,
#'     dat = dat,
#'     id.0 = "G")
#'
#' # Visualize the network
#' plot(network_result$graph)
#' }
#'
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njau.edu.cn}
CNPS.network = function(
    ps = ps,
    dat = dat,
    # path0 = path0,
    id.0 = "C"
){
  # path = paste(path0,"/",id.0,sep = "")
  # fs::dir_create(path)
  id = dat %>%filter(Group == id.0) %>%.$K.number
  tax = dat %>%filter(Group == id.0) %>% distinct(K.number, .keep_all = TRUE)
  head(tax)
  row.names(tax) = tax$K.number

  tax$symbol %>% table()

  ps.t = ps %>%
    scale_micro(method = "rela") %>%
    subset_taxa.wt2("OTU",id)
  tax_table(ps.t) = as.matrix(tax)

  ps.t2 = ps.t %>% tax_glom_wt(ranks = "symbol")

  # ps.t2 = ps.t %>% tax_glom_wt(ranks = "symbol")
  map =sample_data(ps.t2)
  head(map)
  # colnames(map)[colnames(map) == id2[i]] = "Group"
  # gru2 = id2[i]
  # sample_data(ps.t2) = map

  # ps.t3 = ps.t2 %>% subset_samples(zone %in% c(gru1)) %>%
  #   filter_taxa(function(x) sum(x ) > 0, TRUE)
  map.t = ps.t2 %>% sample_data() %>% as.data.frame()


  gru3 = map.t$Group %>% unique() %>% as.character()
  n = 1
  for (n in 1:length(gru3)) {
    ps.t4 = ps.t2 %>% subset_samples.wt("Group",gru3[n]) %>%
      filter_taxa(function(x) sum(x ) > 0, TRUE)
    # ?cor_Big_micro
    result = cor_Big_micro(ps = ps.t4,
                           N = 0,
                           r.threshold=0.6,
                           p.threshold=0.05,
                           method = "pearson"
    )

    #--提取相关矩阵
    cor = result[[1]]
    dim(cor)
    igraph = make_igraph(cor)
    library(igraph)
    # dat = net_properties.4(igraph,n.hub = F)
    head(dat,n = 16)
    # write.csv(dat,paste(path,"/network",id.0,"_",gru3[n],
    #                     "node_properties.csv",sep = ""),quote = F)


    netClu = data.frame(ID = row.names(cor),group =rep(1,length(row.names(cor)))[1:length(row.names(cor))] )
    netClu$group = as.factor(netClu$group)
    result2 = PolygonClusterG (cor = cor,nodeGroup =netClu )
    node = result2[[1]]
    edge = edgeBuild(cor = cor,node = node)
    head(edge)
    if (dim(edge)[1] != 0) {
      edge2 = edge %>% filter(weight != 0)
      tem = c(edge2$OTU_2,edge2$OTU_1) %>% unique()
      cor2 = cor[tem,tem]

      netClu = data.frame(ID = row.names(cor2),group =rep(1,length(row.names(cor2)))[1:length(row.names(cor2))] )
      netClu$group = as.factor(netClu$group)
      result2 = PolygonClusterG (cor = cor2,nodeGroup =netClu )
      node = result2[[1]]
      otu = ps.t4 %>% vegan_otu() %>% t() %>%
        as.data.frame()
      head(otu)
      otu$mean = rowSums(otu)
      otu = otu %>% rownames_to_column("ID")
      nodes = node %>% left_join(otu,by = c("elements"="ID"))
      head(nodes)
      #-----计算边
      edge = edgeBuild(cor = cor2,node = node)
      head(edge)
      edge$abs = abs(edge$weight)

      ### 出图
      library(ggrepel)
      library(ggnewscale)
      p1 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = abs),
                                    data = edge, size = 1) +
        # ggnewscale::new_scale_s() +
        geom_point(aes(X1, X2,size = mean),pch = 21, data = nodes,fill = "#984EA3") +
        scale_color_gradientn(colours =c("white","grey10"))+
        scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
        # labs( title = paste(layout,"network",sep = "_"))+
        ggrepel::geom_text_repel(aes(X1, X2,label=elements),size=4, data = nodes)+
        # discard default grid + titles in ggplot2
        theme(panel.background = element_blank()) +
        # theme(legend.position = "none") +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
        theme(legend.background = element_rect(colour = NA)) +
        theme(panel.background = element_rect(fill = "white",  colour = NA)) +
        theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
      p1
      # write.csv(nodes,paste(path,"/network",id.0,"_",gru3[n],"nodes.csv",sep = ""),quote = F)
      # write.csv(edge,paste(path,"/network",id.0,"_",gru3[n],"nedge.csv",sep = ""),quote = F)
      # ggsave(paste(path,"/network",id.0,"_",gru3[n],".pdf",sep = ""),p1,width = 6,height = 5)
      #
    }
  }
  return(list(p1,nodes,edge))
}
