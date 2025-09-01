#' @title Generate Taxonomic Sankey Diagram for Microbial Data
#'
#' @description
#' This function creates Sankey diagrams to visualize microbial taxonomic transitions from higher to lower taxonomic levels (e.g., Kingdom to Genus).
#' It generates interactive Sankey diagrams for each group in the `phyloseq` object based on the specified taxonomic rank and top OTUs.
#'
#' @param ps A `phyloseq` object containing microbiome data.
#' @param rank A numeric value specifying the taxonomic rank for aggregation. Default is `6` (e.g., Genus level).
#' @param Top An integer specifying the number of top OTUs to include based on abundance. Default is `50`.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{SankeyDiagram}{An interactive Sankey diagram generated with `networkD3::sankeyNetwork`.}
#'   \item{Data}{A data frame containing the Sankey diagram links and node data.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Extracts OTU and taxonomic data for each group in the `phyloseq` object.
#'   \item Aggregates taxonomic information up to the specified rank and filters the top `Top` OTUs by abundance.
#'   \item Creates source-target relationships for taxonomic transitions between levels (e.g., Kingdom → Phylum, Phylum → Class).
#'   \item Calculates mean abundance values for each taxonomic level within groups.
#'   \item Generates an interactive Sankey diagram visualizing taxonomic transitions and their relative abundances.
#' }
#'
#' The Sankey diagram is interactive and allows users to explore taxonomic transitions and their abundances between levels.
#'
#' @examples
#' \dontrun{
#' res = sankey.micro(ps = ps.16s, rank = 6, Top = 50)
#' p22 = res$SankeyDiagram
#' dat = res$Data
#' }
#'
#' @author Contact: Tao Wen \email{2018203048@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njau.edu.cn}
#'
#' @export

sankey.micro = function(ps = ps,
                        rank = 6,# 参数目前不可修改
                        Top = 50

){

id.g = sample_data(ps)$Group %>% unique() %>% as.character()
  print(id.g)
  j = 1
 for (j in 1:length(id.g)) {

   otu = ps %>% vegan_otu() %>% t() %>%
     as.data.frame()
   map = sample_data(ps)
   map$ID = row.names(map)
   sample_data(ps) = map
   otu = otu[,map$ID[as.character(map$Group) == id.g[j]]]
   ps.t = ps

   otu_table(ps.t) = otu_table(as.matrix(otu),taxa_are_rows = TRUE)
   print("1")

   tax = ps.t %>%
     subset_samples.wt("Group", c(id.g[j])) %>%
     ggClusterNet::tax_glom_wt(ranks = rank) %>%
     ggClusterNet::filter_OTU_ps(Top) %>%
     subset_taxa.wt("Genus", c("Unassigned","Unknown"),TRUE) %>%
     ggClusterNet::vegan_tax() %>%
     as.data.frame()
   head(tax)
   dim(tax)
   id2 = c("k","p","c","o","f","g")
   dat = NULL
   for (i in 1:5) {
     dat <- tax[,c(i,i+1)] %>% distinct(.keep_all = TRUE)
     colnames(dat) = c("source","target")
     dat$source = paste(id2[i],dat$source,sep = "_")
     dat$target = paste(id2[i+1],dat$target,sep = "_")
     if (i == 1) {
       dat2 = dat
     }

     dat2 = rbind(dat2,dat)

   }
   dim(dat2)
   # dat2 = dat2 %>% distinct(.keep_all = TRUE)

   head(dat2)
   #print("1")
   otu = ps.t %>%
     subset_samples.wt("Group" , c(id.g[j])) %>%
     ggClusterNet::tax_glom_wt(ranks = 6) %>%
     ggClusterNet::scale_micro() %>%
     ggClusterNet::filter_OTU_ps(Top) %>%
     subset_taxa.wt("Genus" , c("Unassigned","Unknown"),TRUE) %>%
     ggClusterNet::vegan_otu() %>%
     t() %>%
     as.data.frame()

   head(otu)
   otutax = cbind(otu,tax)

   id = rank.names(ps)[1:6]
   dat = NULL

   for (i in 1:6) {
     dat <- otutax %>%
       dplyr::group_by(!!sym(id[i])) %>%
       summarise_if(is.numeric,mean,na.rm = TRUE)

     dat = data.frame(Genus = dat[,1],value =rowSums(dat[,-1]) )
     colnames(dat) = c("target","value")
     dat$target = paste(id2[i],dat$target,sep = "_")
     if (i == 1) {
       dat3 = dat
     }

     dat3 = rbind(dat3,dat)
   }

   dat4 <- dat2 %>% left_join(dat3)
   sankey = dat4

   head( sankey )

   nodes <- data.frame(name = unique(c(as.character(sankey$source),as.character(sankey$target))),stringsAsFactors = FALSE)
   nodes$ID <- 0:(nrow(nodes)-1)
   sankey <- merge(sankey,nodes,by.x = "source",by.y = "name")
   sankey <- merge(sankey,nodes,by.x = "target",by.y = "name")
   colnames(sankey) <- c("X","Y","value","source","target")
   sankey <- subset(sankey,select = c("source","target","value"))
   nodes <- subset(nodes,select = c("name"))

   ColourScal='d3.scaleOrdinal() .range(["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"])'

   sankey$energy_type <- sub(' .*', '', nodes[sankey$source + 1, 'name'])
   library(networkD3)
   p <- sankeyNetwork(Links = sankey, Nodes = nodes,
                      Source = "source",Target = "target",Value = "value",
                      NodeID = "name",
                      sinksRight=FALSE,
                      LinkGroup = 'energy_type',
                      colourScale= ColourScal,
                      # nodeWidth=40,
                      # fontSize=13
                      # nodePadding=20
   )
   # i = 1
   dat = sankey

   # subsnapath = paste(snapath,"/",id.g[j],sep = "")
   # dir.create(subsnapath)
   # FileName <-paste(subsnapath,"/",id.g[j],"sankey1.csv", sep = "")
   # write.csv(dat,FileName,sep = "")

   # saveNetwork(p,paste(subsnapath,"/sankey1.html", sep = ""))
   # library(webshot)
   # # webshot::install_phantomjs()
   # # webshot(paste(snapath,"/sankey1.html", sep = "") ,paste(snapath,"/sankey1.png", sep = ""))
   # webshot(paste(subsnapath,"/sankey1.html", sep = "") , paste(subsnapath,"/sankey1.pdf", sep = ""))
   #
   return(list(p,dat))
 }

}
