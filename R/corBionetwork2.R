corBionetwork2 = function (otu = NULL, tax = NULL, map = NULL, ps = NULL, lab = NULL,
          N = 0, r.threshold = 0.6, p.threshold = 0.05, label = FALSE,
          group = "Group", env = NULL, envGroup = NULL, method = "spearman",
          layout = "fruchtermanreingold", path = "./", fill = "Phylum",
          size = "igraph.degree", scale = TRUE, bio = TRUE, zipi = FALSE,
          step = 100, width = 20, height = 20, big = TRUE, select_layout = TRUE,
          layout_net = "model_maptree", clu_method = "cluster_fast_greedy",
          minsize = 4, maxsize = 14)
{
 # dir.create(path)
  ps = ggClusterNet::inputMicro(otu, tax, map, tree, ps, group = group)
  if (scale) {
    ps = ps %>% ggClusterNet::scale_micro()
  }
  ps_all = ps %>% ggClusterNet::filter_OTU_ps(N)
  mapping = as.data.frame(phyloseq::sample_data(ps))
  mapping$ID = row.names(mapping)
  sample_data(ps) = mapping
  y = matrix(1, nrow = 14, ncol = length(unique(mapping$Group)))
  layouts = as.character(unique(mapping$Group))
  aa = 1
  plots = list()
  layout = layouts[1]
  layout
  for (layout in layouts) {
    print(layout)
    map <- as.data.frame(phyloseq::sample_data(ps))
    mapsub <- map[map$Group == layout, ]
    ps_sub <- ps
    sample_data(ps_sub) <- mapsub
    ps_sub = phyloseq::filter_taxa(ps_sub, function(x) sum(x) >  0, TRUE)
    if (!is.null(env)) {
      colnames(env)[1] = "ID"
      env_sub <- env[match(mapsub$ID, env$ID), ]
      head(env_sub)
    }
    if (bio) {
      if (!is.null(env)) {
        if (big) {
          result <- corBiostripeBig(data = env_sub, group = envGroup,
                                    ps = ps_sub, r.threshold = r.threshold, p.threshold = p.threshold,
                                    method = method)
        }
        else {
          result <- corBiostripe(data = env_sub, group = envGroup,
                                 ps = ps_sub, r.threshold = r.threshold, p.threshold = p.threshold,
                                 method = method)
        }
        occor.r = result[[1]]
        tax = as.data.frame((ggClusterNet::vegan_tax(ps_sub)))
        if (length(tax$filed) != 0) {
          group2 <- data.frame(SampleID = row.names(tax),
                               Group = tax$filed)
        }
        else {
          group2 <- data.frame(SampleID = row.names(tax),
                               Group = "OTU")
        }
        colnames(envGroup) <- c("SampleID", "Group")
        netClu = rbind(envGroup, group2)
        colnames(netClu) <- c("ID", "group")
      }
      else {
        if (big) {
          result <- corBiostripeBig(ps = ps_sub, r.threshold = r.threshold,
                                    p.threshold = p.threshold, method = method)
        }
        else {
          result <- corBiostripe(ps = ps_sub, r.threshold = r.threshold,
                                 p.threshold = p.threshold, method = method)
        }
        occor.r = result[[1]]
        tax = as.data.frame((ggClusterNet::vegan_tax(ps_sub)))
        if (length(tax$filed) != 0) {
          group2 <- data.frame(SampleID = row.names(tax),
                               Group = tax$filed)
        }
        else {
          group2 <- data.frame(SampleID = row.names(tax),
                               Group = "OTU")
        }
        netClu = group2
        colnames(netClu) <- c("ID", "group")
      }
    }
    result4 = ggClusterNet::nodeEdge(cor = occor.r)
    igraph = igraph::graph_from_data_frame(result4[[1]],
                                           directed = FALSE, vertices = result4[[2]])
    if (zipi) {
      print("zipi_start")
      res = ZiPiPlot(igraph = igraph, method = "cluster_fast_greedy")
      p <- res[[1]]
      ggsave(paste(path, "/", layout, "_ZiPi.pdf", sep = ""),
             p)
      ZiPi <- res[[2]]
      write.csv(ZiPi, paste(path, "/", layout, "ZiPi.csv",
                            sep = ""), row.names = FALSE)
    }
    netClu$group = as.factor(netClu$group)
    if (select_layout) {
      node = NULL
      nrow(occor.r)

      node = culculate_node_axis(cor.matrix = occor.r,
                                 layout = layout_net, seed = 1, group = NULL,
                                 model = FALSE, method = clu_method)
    }
    else if (select_layout) {
      result2 <- model_Gephi.2(cor = cor, method = clu_method,
                               seed = 12)
      node = result2[[1]]
    }
    tax = ps_sub %>% ggClusterNet::vegan_tax() %>% as.data.frame()
    nodesub1 <- merge(node, tax, by = "row.names", all = T)
    row.names(nodesub1) = nodesub1$Row.names
    nodesub1$Row.names = NULL
    plotcord <- nodesub1 %>% dplyr::inner_join(netClu, by = c(elements = "ID"))
    edges = ggClusterNet::edgeBuild(cor = occor.r, node = node)
    head(edges)
    edge_Gephi = data.frame(source = edges$OTU_1, target = edges$OTU_2,
                            correlation = edges$weight, direct = "undirected",
                            cor = edges$cor)
    node_Gephi = data.frame(ID = plotcord$elements, plotcord[4:dim(plotcord)[2]],
                            Label = plotcord$elements)
    write.csv(edge_Gephi, paste(path, "/", layout, "_Gephi_edge.csv",
                                sep = ""), row.names = FALSE)
    write.csv(node_Gephi, paste(path, "/", layout, "_Gephi_node.csv",
                                sep = ""), row.names = FALSE)
    nodepro = ggClusterNet::node_properties(igraph)
    write.csv(nodepro, paste(path, "/", layout, "_node_properties.csv",
                             sep = ""), row.names = FALSE)
    row.names(plotcord) = plotcord$elements
    nodeG = merge(plotcord, nodepro, by = "row.names", all.x = T)
    row.names(nodeG) = nodeG$Row.names
    nodeG$Row.names = NULL
    numna = (dim(nodeG)[2] - 3):dim(nodeG)[2]
    nodeG[, numna][is.na(nodeG[, numna])] = 0
    head(nodeG)
    p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2,
                                      yend = Y2, color = as.factor(cor)), data = edges,
                                  size = 0.3, alpha = 0.5) + geom_point(aes(x = X1,
                                                                            y = X2, size = igraph.degree, fill = group), pch = 21,
                                                                        data = nodeG) + scale_colour_brewer(palette = "Set1") +
      scale_size(range = c(minsize, maxsize)) + scale_x_continuous(breaks = NULL) +
      scale_y_continuous(breaks = NULL) + labs(title = paste(layout,
                                                             "network", sep = "_")) + theme_void()
    p0
    head(nodeG)
    tem = nodeG %>% dplyr::filter(elements %in% lab[[1]])
    if (label == TRUE) {
      if (!is.null(lab)) {
        p1 <- p0 + ggrepel::geom_text_repel(aes(X1, X2,
                                                label = elements), size = 4, data = tem)
      }
      else {
        p1 <- p0 + ggrepel::geom_text_repel(aes(X1, X2,
                                                label = elements), size = 4, data = nodeG)
      }
      plotname = paste(path, "/network_lab_", layout, ".pdf",
                       sep = "")
      ggsave(plotname, p1, width = width, height = height)
    }
    p1
    plotname = paste(path, "/network", layout, ".pdf", sep = "")
    ggsave(plotname, p0, width = width, height = height)
    print("1")
    plots[[aa]] = p0
    rand.g <- erdos.renyi.game(length(V(igraph)), length(E(igraph)),
                               type = c("gnm"))
    data1 = data.frame(network = degree_distribution(igraph,
                                                     cumulative = FALSE), group = "Erdős–Rényi network",
                       ID = c(1:length(degree_distribution(igraph, cumulative = FALSE))))
    data2 = data.frame(network = degree_distribution(rand.g,
                                                     cumulative = FALSE), group = "network", ID = c(1:length(degree_distribution(rand.g,
                                                                                                                                 cumulative = FALSE))))
    data = rbind(data1, data2)
    p1 <- ggplot(data) + geom_point(aes(x = ID, y = network,
                                        group = group, fill = group), pch = 21, size = 2) +
      geom_smooth(aes(x = ID, y = network, group = group,
                      color = group)) + theme_bw()
    plotname = paste(path, "/Power_law_distribution_", layout,
                     ".pdf", sep = "")
    ggsave(plotname, p1, width = width, height = height)
    rand.g.netpro_result <- c()
    for (i in 1:step) {
      rand.g <- erdos.renyi.game(length(V(igraph)), length(E(igraph)),
                                 type = c("gnm"))
      tem_netpro_result <- ggClusterNet::net_properties(rand.g)
      rand.g.netpro_result <- cbind(rand.g.netpro_result,
                                    tem_netpro_result)
    }
    print("2")
    result_summary <- cbind(rowMeans(rand.g.netpro_result),
                            apply(rand.g.netpro_result, 1, sd))
    colnames(result_summary) <- c("Means", "SD")
    igraph.weight <- igraph::E(igraph)$weight
    E(igraph)$weight <- NA
    igraph <- igraph::remove.edge.attribute(igraph, "weight")
    netpro_result <- ggClusterNet::net_properties(igraph)
    colnames(netpro_result) <- layout
    sum_net = cbind(netpro_result, result_summary)
    write.csv(sum_net, paste(path, "/", layout, "_net_VS_erdos_properties.csv",
                             sep = ""), row.names = TRUE)
    print("3")
    y = as.data.frame(y)
    colnames(y) = layouts
    y[layout] = netpro_result[, 1]
    row.names(y) = row.names(netpro_result)
    aa = aa + 1
  }
  plotname = paste(path, "/network_all.pdf", sep = "")
  p = ggpubr::ggarrange(plotlist = plots, common.legend = TRUE,
                        legend = "right")
  return(list(p, y, edges, nodeG, occor.r))
}
