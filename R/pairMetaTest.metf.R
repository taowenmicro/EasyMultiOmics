#' @title Beta diversity statistics by adonis/anosim/MRPP in pair
#' @description Input phyloseq object, test method and distance type.
#' @param otu Metagenome functional composition table.
#' @param map Sample metadata;
#' @param ps A phyloseq format file used as an alternative for the input containing metagenome functional composition table, tax, and sample metadata.
#' @param tree Tree/nwk file;
#' @param dist Distance type, including "unifrac" "wunifrac" "dpcoa" "jsd" "manhattan" "euclidean"   "canberra" "bray" "kulczynski"  "jaccard" "gower" "altGower" "morisita" "horn" "mountford"  "raup" "binomial"  "chao"  "cao" "w"  "-1"  "c" "wb"  "r"   "I"  "e" "t" "me"   "j"  "sor"  "m"   "-2"  "co";
#' @param group Column name for groupID in map table(sample metadata).
#' @param method DCA, CCA, RDA, NMDS, MDS, PCoA, PCA, LDA;
#' @param pvalue.cutoff Pvalue threshold, default in 0.05.
#' @param Micromet Statistics default by adonis, alternative anosim or MRPP;
#' @details
#' By default, input phyloseq object include metadata and metagenome functional composition table
#' The available diversity indices include the following:
#' \itemize{
#' \item{most used indices: bray unifrac wunifrac}
#' \item{other used indices: dpcoa jsd manhattan euclidean canberra kulczynski jaccard gower altGower morisita horn mountford raup binomial chao cao w -1 c wb r I e t me j sor m -2 co}
#' }
#' @return stat table
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' ps =ps.kegg %>% filter_OTU_ps(Top = 1000)
#' dat2 = pairMicroTest.metf(ps = ps, Micromet = "MRPP", dist = "bray")
#' dat2
#' @export
pairMetaTest.metf = function(ps = ps, Micromet = "adonis", dist = "bray"){

  ps1_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) );ps1_rela

  map = as.data.frame(phyloseq::sample_data(ps1_rela))

  otu = ps1_rela %>% ggClusterNet::vegan_otu() %>%
    as.data.frame()
  map$Group = as.factor(map$Group)
  aa = levels(map$Group)
  # aa
  aaa = combn(aa,2)
  aaa
  dim(aaa)[2]

  # 构建三个空列
  ID = rep("a",dim(aaa)[2])
  R = rep("a",dim(aaa)[2])
  P = rep("a",dim(aaa)[2])

  for (i in 1:dim(aaa)[2]) {
    print(i)
    Desep_group = aaa[,i]
    map = as.data.frame(phyloseq::sample_data(ps1_rela))
    # head(map)
    map$ID = row.names(map)

    # maps<- dplyr::filter(map,Group %in% Desep_group)
    maps <- map[map$Group %in%Desep_group,]
    row.names(maps) = maps$ID
    ps_sub = ps1_rela
    phyloseq::sample_data( ps_sub ) = maps
    ps_sub = phyloseq::filter_taxa(ps_sub, function(x) sum(x ) > 0 , TRUE);ps_sub
    map = as.data.frame(phyloseq::sample_data(ps_sub))
    gg =map$Group
    unif <- phyloseq::distance(ps_sub, method=dist)

    if (Micromet == "MRPP") {
      mrpp = vegan::mrpp(unif, map$Group)
      as1 = round(mrpp$delta,3)
      R2 <- paste("MRPP.delta ",as1, sep = "")
      # R[i] = R2
      R2
      p_v = paste("p: ",round(mrpp$Pvalue,3), sep = "")
      p_v
    }

    if (Micromet == "anosim") {
      dat.ano = vegan::anosim(unif, map$Group)
      a = round(dat.ano$statistic,3)
      R2 <- paste("ANOSIM.r ",a, sep = "")
      R[i] = R2
      p_v = paste("p: ",round(dat.ano$signif,3), sep = "")
      P[i] = p_v
    }
    if (Micromet == "adonis") {

      # map_df <- shiny::reactive({map})
      str(otu)
      str(dune.env)
      ado = vegan:: adonis2( unif~ map$Group,method = "bray", by = NULL)
      R2 = paste("Adonis:R ",round(ado$R2[1],3), sep = "")
      p_v = paste("p: ",ado$`Pr(>F)`[1], sep = "")
      R[i] = R2
      P[i] = p_v
      # title = paste(R2," ",p_v, sep = "")
      # title
      # print(i)
      # print(R)
    }
    ID[i] = paste(Desep_group[1],Desep_group[2],sep = "_VS_")
    P[i] = p_v
    R[i] = R2
  }

  result = data.frame(ID = ID,stat = R,p = P)


  return(result)
}
