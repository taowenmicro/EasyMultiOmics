
mental.omics= function(ps01= ps01,ps02= ps02,method="spearman" ){
dist.01 = ps01 %>%
  scale_micro() %>%
  # tax_glom_wt(j) %>%
  # subset_taxa.wt("OTU",id.micro) %>%
  otu_table() %>% t() %>%
  vegan::vegdist(method="bray") %>%
  as.matrix()




dist.02 = ps02 %>%
  scale_micro() %>%
  # subset_taxa.wt("OTU",id.ms) %>%
  otu_table() %>%
  t() %>%
  vegan::vegdist(method="bray") %>%
  as.matrix()

# method =  "spearman"
# method =  "pearson"

dist.list = list(Bac =dist.01,
                 Fun = dist.02
)

sample_data(ps01)

gru = c("da1","da2")
id = combn(unique(gru),2)
names(dist.list) = gru


R_mantel = c()
p_mantel = c()
name = c()
R_pro <- c()
p_pro <- c()
plots = list()
i = 1
for (i in 1:dim(id)[2]) {


  dist1 = dist.list[[id[1,i]]]
  dist2 <- dist.list[[id[2,i]]]
  mt <- vegan::mantel(dist1,dist2,method = method)
  R_mantel[i] = mt$statistic
  p_mantel[i] = mt$signif

  name[i] = paste(id[1,i],"_VS_",id[2,i],sep = "")
  #--p

  mds.s <- vegan::monoMDS(dist1)
  mds.r <- vegan::monoMDS(dist2)
  pro.s.r <- vegan::protest(mds.s,mds.r)

  R_pro[i] <- pro.s.r$ss
  p_pro[i] <- pro.s.r$signif

}
dat = data.frame(name,R_mantel,p_mantel,R_pro,p_pro )
head(dat)
return(dat)
}
