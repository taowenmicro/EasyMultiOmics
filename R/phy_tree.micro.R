# library(ggtreeExtra)
library(ggtree)
# library(treeio)
# library(ggstar)
#创建新的scale，多个fill或者color时
library(ggnewscale)
# #-读取进化树
# library(patchwork)
# library(ggClusterNet)
# library(phyloseq)
# library(tidyverse)
#--微生物组进化树功能
#
# barpath = paste(otupath,"/phy_tree_micro/",sep = "")
# dir.create(barpath)
#


remove_rankID = function(taxtab){
  taxtab$Kingdom = gsub("d__","",taxtab$Kingdom)
  taxtab$Kingdom = gsub("k__","",taxtab$Kingdom)
  taxtab$Phylum = gsub("p__","",taxtab$Phylum)
  taxtab$Class = gsub("c__","",taxtab$Class)
  taxtab$Order = gsub("o__","",taxtab$Order)
  taxtab$Family = gsub("f__","",taxtab$Family)
  taxtab$Genus = gsub("g__","",taxtab$Genus)
  taxtab$Species = gsub("s__","",taxtab$Species)
  return(taxtab)
}

# detach("package:treeio")
# detach("package:vegan")

#
# result <- phy_tree_micro (ps = ps,Top = 100)
#
# p0 = result[[1]]
# p1 = result[[2]]
# p2 = result[[3]]
# p3 = result[[4]]
# p4 = result[[5]]
#
#
#
# FileName <- paste(barpath,Top_micro,"phy_tree_micro1", ".pdf", sep = "")
# ggsave(FileName, p0, width = 6, height = 6)
#
# FileName <- paste(barpath,Top_micro,"phy_tree_micro2", ".pdf", sep = "")
# ggsave(FileName, p1, width = 7, height = 7)
#
# FileName <- paste(barpath,Top_micro,"phy_tree_micro3", ".pdf", sep = "")
# ggsave(FileName, p2, width = 7, height = 7)
#
# FileName <- paste(barpath,Top_micro,"phy_tree_micro4", ".pdf", sep = "")
# ggsave(FileName, p3, width = 12, height = 12)
#
# FileName <- paste(barpath,Top_micro,"phy_tree_micro5", ".pdf", sep = "")
# ggsave(FileName, p4, width = 15, height = 15)

phy_tree.micro2 <- function(
  ps = ps,
  Top = 100
){
  tax = ps %>% vegan_tax() %>%
    as.data.frame()
  head(tax)
  tax = remove_rankID(tax) %>%as.matrix()
  tax[is.na(tax)] = "Unknown"
  tax[tax == " "] = "Unknown"
  tax_table(ps) = as.matrix(tax)

  alltax = ps %>%
    ggClusterNet::filter_OTU_ps(Top) %>%
    ggClusterNet::vegan_tax() %>%
    as.data.frame()
  alltax$OTU = row.names(alltax)

  head(alltax)
  # alltax$Kingdom = sapply(strsplit(alltax$Kingdom, ":"), `[`, 2)
  # alltax$Phylum = sapply(strsplit(alltax$Phylum, ":"), `[`, 2)
  # alltax$Class = sapply(strsplit(alltax$Class, ":"), `[`, 2)
  # alltax$Order = sapply(strsplit(alltax$Order, ":"), `[`, 2)
  # alltax$Family = sapply(strsplit(alltax$Family, ":"), `[`, 2)
  # alltax$Genus = sapply(strsplit(alltax$Genus, ":"), `[`, 2)
  # alltax$Species = sapply(strsplit(alltax$Species, ":"), `[`, 2)
  # alltax[is.na(alltax)] = "Unknown"

  trda <- MicrobiotaProcess::convert_to_treedata(alltax)
  p0 <- ggtree::ggtree(trda, layout="inward_circular", size=0.2, xlim=c(30,NA))

  # p0$data
  tax = ps %>%
    ggClusterNet::filter_OTU_ps(Top) %>%
    ggClusterNet::vegan_tax() %>%
    as.data.frame()

  tippoint = data.frame(OTU = row.names(tax),Taxa = tax$Phylum,Level = "Phylum")

  tippoint$OTU = paste("",tippoint$OTU,sep = "")
  tippoint$names <- gsub("","",tippoint$OTU)


  p0_1 <- ggtree::ggtree(trda, layout="circular", size=0.2, xlim=c(30,NA)) %<+% tippoint
  p0_1
  p0 <- p0 %<+% tippoint
  p0
  a <- tippoint$Taxa %>% unique() %>% length()

  b = rep(18,a)
  names(b) = tippoint$Taxa %>% unique()


  p1 <- p0 +
    geom_tippoint(
      mapping=aes(
        color=Taxa,
        shape=Level
      ),
      size=1,
      alpha=0.8
    ) +
    scale_color_manual(values=colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(a),
                       guide=guide_legend(
                         keywidth=0.5,
                         keyheight=0.5,
                         order=2,
                         override.aes=list(shape= b,
                                           size=2
                         ),
                         na.translate=TRUE
                       )
    ) +
    scale_shape_manual(values=c("Phylum"=20, "Class"=18), guide="none" )


  p1_1 <- p0_1 +
    geom_tippoint(
      mapping=aes(
        color=Taxa,
        shape=Level
      ),
      size=1,
      alpha=0.8
    ) +
    scale_color_manual(values=colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(a),
                       guide=guide_legend(
                         keywidth=0.5,
                         keyheight=0.5,
                         order=2,
                         override.aes=list(shape= b,
                                           size=2
                         ),
                         na.translate=TRUE
                       )
    ) +
    scale_shape_manual(values=c("Phylum"=20, "Class"=18), guide="none" )


  #------对OTU之间的关系进行连线
  otu = ps %>%
    ggClusterNet::filter_OTU_ps(Top) %>%
    ggClusterNet::vegan_otu() %>%
    t() %>%
    as.data.frame()
  head(otu)
  pssub = ps %>%
    ggClusterNet::filter_OTU_ps(Top)
  result = ggClusterNet::corMicro (ps = pssub ,N = 0,r.threshold=0.8,p.threshold=0.05,method = "pearson")
  #--提取相关矩阵
  cor = result[[1]]
  diag(cor) = 0

  # library(tidyfst)
  linktab = tidyfst::mat_df(cor) %>%
    dplyr::filter(value != 0) %>%
    dplyr::mutate(direct = ifelse(value> 0, "positive", "nagetive"))
  colnames(linktab) = c("Inhibitor","Sensitive","Interaction","direct")
  head(linktab)
  linktab$Inhibitor = paste("",linktab$Inhibitor,sep = "")
  linktab$Sensitive = paste("",linktab$Sensitive,sep = "")
  head(linktab)

  # p1$data
  p2 <- p1 +
    ggnewscale::new_scale_color() +
    ggtree::geom_taxalink(
      data=linktab,
      mapping=aes(
        taxa1=Inhibitor,
        taxa2=Sensitive,
        color=direct
      ),
      alpha=0.6,
      offset=0.1,
      size=0.15,
      ncp=10,
      hratio=1,
      arrow=grid::arrow(length = unit(0.005, "npc"))
    ) +
    scale_colour_manual(values=c("chocolate2", "#3690C0", "#009E73"),
                        guide=guide_legend(
                          keywidth=0.8, keyheight=0.5,
                          order=1, override.aes=list(alpha=1, size=0.5)
                        )
    )



  otu$id = row.names(otu)
  ringdat = otu %>% tidyfst::longer_dt(id)
  ringdat$value = log2(ringdat$value+1)
  ringdat$id = paste("",ringdat$id,sep = "")
  head(ringdat)
  num <- ringdat$name %>% unique() %>% length()
  p3 <- p2 +
    geom_fruit(
      data=ringdat,
      geom=geom_star,
      mapping=aes(
        y=id,
        x=name,
        size=value,
        fill= name
      ),
      starshape = 13,
      starstroke = 0,
      offset=-0.9,
      pwidth=0.8,
      grid.params=list(linetype=3)
    ) +
    scale_size_continuous(range=c(0, 2),
                          limits=c(sort(ringdat$value)[2], max(ringdat$value)),
                          breaks=c(1, 2, 3),
                          name=bquote(paste(Log[2],"(",.("Count+1"), ")")),
                          guide=guide_legend(keywidth = 0.4, keyheight = 0.4, order=4,
                                             override.aes = list(starstroke=0.3))
    ) +
    scale_fill_manual(
      values=colorRampPalette(RColorBrewer::brewer.pal(12,"Spectral"))(num),
      guide=guide_legend(keywidth = 0.4, keyheight = 0.4, order=3)
    )


  p4 <- p3 +
    geom_tiplab(
      mapping=aes(
        label=names
      ),
      align=TRUE,
      size=2,
      linetype=NA,
      offset=16
    )


  map = phyloseq::sample_data(pssub)

  dat <- pssub %>%
    ggClusterNet::vegan_otu() %>%
    as.data.frame()
  bartab = cbind(dat,data.frame(row.names = row.names(map),ID = row.names(map),Group = map$Group )) %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise_if(is.numeric,mean) %>%
    as.data.frame()%>%
    tidyfst::longer_dt(Group)
  head(bartab)
  bartab$name = paste("",bartab$name,sep = "")
  num = bartab$Group %>% unique() %>%
    length()
  p5 <- p4 +
    ggnewscale::new_scale_fill() +
    geom_fruit(
      data=bartab,
      geom=geom_bar,
      mapping=aes(
        x=value,
        y=name,
        fill= Group
      ),
      stat="identity",
      orientation="y",
      offset=0.48,
      pwidth=1.5,
      axis.params=list(
        axis = "x",
        text.angle = -45,
        hjust = 0,
        vjust = 0.5,
        nbreak = 4
      )
    ) +
    scale_fill_manual(
      name = "Number of interactions",
      values=colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(num),
      guide=guide_legend(keywidth=0.5,keyheight=0.5,order=5)
    ) +
    theme(
      legend.background=element_rect(fill=NA),
      legend.title=element_text(size=6.5),
      legend.text=element_text(size=5),
      legend.spacing.y = unit(0.02, "cm"),
      legend.margin=ggplot2::margin(0.1, 0.9, 0.1,-0.9, unit="cm"),
      legend.box.margin=ggplot2::margin(0.1, 0.9, 0.1, -0.9, unit="cm"),
      plot.margin = unit(c(-1.2, -1.2, -1.2, 0.1),"cm")
    )


  p5_1 <- p1_1 +
    new_scale_fill() +
    geom_fruit(
      data=bartab,
      geom=geom_bar,
      mapping=aes(
        x=value,
        y=name,
        fill= Group
      ),
      stat="identity",
      orientation="y",
      offset=0.1,
      pwidth=1.5,
      axis.params=list(
        axis = "x",
        text.angle = -45,
        hjust = 0,
        vjust = 0.5,
        nbreak = 4
      )
    ) +
    scale_fill_manual(
      name = "Number of interactions",
      values=colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(num),
      guide=guide_legend(keywidth=0.5,keyheight=0.5,order=5)
    ) +
    theme(
      legend.background=element_rect(fill=NA),
      legend.title=element_text(size=6.5),
      legend.text=element_text(size=5),
      legend.spacing.y = unit(0.02, "cm"),
      legend.margin=ggplot2::margin(0.1, 0.9, 0.1,-0.9, unit="cm"),
      legend.box.margin=ggplot2::margin(0.1, 0.9, 0.1, -0.9, unit="cm"),
      plot.margin = unit(c(-1.2, -1.2, -1.2, 0.1),"cm")
    )


  return(list(p0,p1,p2,p3,p4,p5,p1_1,p5_1))
}


phy_tree.micro <- function(
  ps = ps,
  Top = 100
){
  tax = ps %>% vegan_tax() %>%
    as.data.frame()
  head(tax)
  tax = remove_rankID(tax) %>%as.matrix()
  tax[is.na(tax)] = "Unknown"
  tax[tax == " "] = "Unknown"
  tax_table(ps) = as.matrix(tax)
  alltax = ps %>%
    ggClusterNet::filter_OTU_ps(Top) %>%
    ggClusterNet::vegan_tax() %>%
    as.data.frame()
  alltax$OTU = row.names(alltax)
  head(alltax)
  # alltax$Kingdom = sapply(strsplit(alltax$Kingdom, ":"), `[`, 2)
  # alltax$Phylum = sapply(strsplit(alltax$Phylum, ":"), `[`, 2)
  # alltax$Class = sapply(strsplit(alltax$Class, ":"), `[`, 2)
  # alltax$Order = sapply(strsplit(alltax$Order, ":"), `[`, 2)
  # alltax$Family = sapply(strsplit(alltax$Family, ":"), `[`, 2)
  # alltax$Genus = sapply(strsplit(alltax$Genus, ":"), `[`, 2)
  # alltax$Species = sapply(strsplit(alltax$Species, ":"), `[`, 2)
  # alltax[is.na(alltax)] = "Unknown"


  trda <-  MicrobiotaProcess::convert_to_treedata(alltax)
  p0 <- ggtree(trda, layout="inward_circular", size=0.2, xlim=c(30,NA))
  # p0$data
  tax = ps %>%
    ggClusterNet::filter_OTU_ps(Top) %>%
    ggClusterNet::vegan_tax() %>%
    as.data.frame()

  tippoint = data.frame(OTU = row.names(tax),Taxa = tax$Phylum,Level = "Phylum")

  tippoint$OTU = paste("st__",tippoint$OTU,sep = "")
  tippoint$names <- gsub("st__","",tippoint$OTU)

  head(tippoint)

  p0_1 <- ggtree(trda, layout="circular", size=0.2, xlim=c(30,NA)) %<+% tippoint
  p0_1
  p0 <- p0 %<+% tippoint
  p0
  a <- tippoint$Taxa %>% unique() %>% length()

  b = rep(18,a)
  names(b) = tippoint$Taxa %>% unique()


  p1 <- p0 +
    geom_tippoint(
      mapping=aes(
        color=Taxa,
        shape=Level
      ),
      size=1,
      alpha=0.8
    ) +
    scale_color_manual(values=colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(a),
                       guide=guide_legend(
                         keywidth=0.5,
                         keyheight=0.5,
                         order=2,
                         override.aes=list(shape= b,
                                           size=2
                         ),
                         na.translate=TRUE
                       )
    ) +
    scale_shape_manual(values=c("Phylum"=20, "Class"=18), guide="none" )

  p0_1
  p1_1 <- p0_1 +
    geom_tippoint(
      mapping=aes(
        color=Taxa,
        shape=Level
      ),
      size=1,
      alpha=0.8
    ) +
    scale_color_manual(values=colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(a),
                       guide=guide_legend(
                         keywidth=0.5,
                         keyheight=0.5,
                         order=2,
                         override.aes=list(shape= b,
                                           size=2
                         ),
                         na.translate=TRUE
                       )
    ) +
    scale_shape_manual(values=c("Phylum"=20, "Class"=18), guide="none" )

  p1_1
  #------对OTU之间的关系进行连线
  otu = ps %>%
    ggClusterNet::filter_OTU_ps(Top) %>%
    ggClusterNet::vegan_otu() %>%
    t() %>%
    as.data.frame()
  head(otu)
  pssub = ps %>%
    ggClusterNet::filter_OTU_ps(Top)
  result = ggClusterNet::corMicro (ps = pssub ,N = 0,r.threshold=0.8,p.threshold=0.05,method = "pearson")
  #--提取相关矩阵
  cor = result[[1]]
  diag(cor) = 0

  # library(tidyfst)
  linktab = tidyfst::mat_df(cor) %>%
    dplyr::filter(value != 0) %>%
    dplyr::mutate(direct = ifelse(value> 0, "positive", "nagetive"))
  colnames(linktab) = c("Inhibitor","Sensitive","Interaction","direct")
  head(linktab)
  linktab$Inhibitor = paste("st__",linktab$Inhibitor,sep = "")
  linktab$Sensitive = paste("st__",linktab$Sensitive,sep = "")
  head(linktab)


  p2 <- p1 +
    ggnewscale::new_scale_color() +
    ggtree::geom_taxalink(
      data=linktab,
      mapping=aes(
        taxa1=Inhibitor,
        taxa2=Sensitive,
        color=direct
      ),
      alpha=0.6,
      offset=0.1,
      size=0.15,
      ncp=10,
      hratio=1,
      arrow=grid::arrow(length = unit(0.005, "npc"))
    ) +
    scale_colour_manual(values=c("chocolate2", "#3690C0", "#009E73"),
                        guide=guide_legend(
                          keywidth=0.8, keyheight=0.5,
                          order=1, override.aes=list(alpha=1, size=0.5)
                        )
    )


p2
  otu$id = row.names(otu)
  ringdat = otu %>% tidyfst::longer_dt(id)
  ringdat$value = log2(ringdat$value+1)
  ringdat$id = paste("st__",ringdat$id,sep = "")
  head(ringdat)
  num <- ringdat$name %>% unique() %>% length()
  p3 <- p2 +
    ggtreeExtra::geom_fruit(
      data=ringdat,
      geom=geom_star,
      mapping=aes(
        y=id,
        x=name,
        size=value,
        fill= name
      ),
      starshape = 13,
      starstroke = 0,
      offset=-0.9,
      pwidth=0.8,
      grid.params=list(linetype=3)
    ) +
    scale_size_continuous(range=c(0, 2),
                          limits=c(sort(ringdat$value)[2], max(ringdat$value)),
                          breaks=c(1, 2, 3),
                          name=bquote(paste(Log[2],"(",.("Count+1"), ")")),
                          guide=guide_legend(keywidth = 0.4, keyheight = 0.4, order=4,
                                             override.aes = list(starstroke=0.3))
    ) +
    scale_fill_manual(
      values=colorRampPalette(RColorBrewer::brewer.pal(12,"Spectral"))(num),
      guide=guide_legend(keywidth = 0.4, keyheight = 0.4, order=3)
    )

  p3
  p4 <- p3 +
    geom_tiplab(
      mapping=aes(
        label=names
      ),
      align=TRUE,
      size=2,
      linetype=NA,
      offset=16
    )


  map = phyloseq::sample_data(pssub)

  dat <- pssub %>%
    ggClusterNet::vegan_otu() %>%
    as.data.frame()
  bartab = cbind(dat,data.frame(row.names = row.names(map),ID = row.names(map),Group = map$Group )) %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise_if(is.numeric,mean) %>%
    as.data.frame()%>%
    tidyfst::longer_dt(Group)

  bartab$name = paste("st__",bartab$name,sep = "")
  num = bartab$Group %>% unique() %>%
    length()
  p4
  p5 <- p4 +
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(
      data=bartab,
      geom=geom_bar,
      mapping=aes(
        x=value,
        y=name,
        fill= Group
      ),
      stat="identity",
      orientation="y",
      offset=0.48,
      pwidth=1.5,
      axis.params=list(
        axis = "x",
        text.angle = -45,
        hjust = 0,
        vjust = 0.5,
        nbreak = 4
      )
    ) +
    scale_fill_manual(
      name = "Number of interactions",
      values=colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(num),
      guide=guide_legend(keywidth=0.5,keyheight=0.5,order=5)
    ) +
    theme(
      legend.background=element_rect(fill=NA),
      legend.title=element_text(size=6.5),
      legend.text=element_text(size=5),
      legend.spacing.y = unit(0.02, "cm"),
      legend.margin=ggplot2::margin(0.1, 0.9, 0.1,-0.9, unit="cm"),
      legend.box.margin=ggplot2::margin(0.1, 0.9, 0.1, -0.9, unit="cm"),
      plot.margin = unit(c(-1.2, -1.2, -1.2, 0.1),"cm")
    )

  p5

  p5_1 <- p1_1 +
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(
      data=bartab,
      geom=geom_bar,
      mapping=aes(
        x=value,
        y=name,
        fill= Group
      ),
      stat="identity",
      orientation="y",
      offset=0.1,
      pwidth=1.5,
      axis.params=list(
        axis = "x",
        text.angle = -45,
        hjust = 0,
        vjust = 0.5,
        nbreak = 4
      )
    ) +
    scale_fill_manual(
      name = "Number of interactions",
      values=colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(num),
      guide=guide_legend(keywidth=0.5,keyheight=0.5,order=5)
    ) +
    theme(
      legend.background=element_rect(fill=NA),
      legend.title=element_text(size=6.5),
      legend.text=element_text(size=5),
      legend.spacing.y = unit(0.02, "cm"),
      legend.margin=ggplot2::margin(0.1, 0.9, 0.1,-0.9, unit="cm"),
      legend.box.margin=ggplot2::margin(0.1, 0.9, 0.1, -0.9, unit="cm"),
      plot.margin = unit(c(-1.2, -1.2, -1.2, 0.1),"cm")
    )
  return(list(p0,p1,p2,p3,p4,p5,p1_1,p5_1,tippoint))
}

