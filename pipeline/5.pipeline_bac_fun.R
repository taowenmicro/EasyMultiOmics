#细菌-真菌联合------
rm(list=ls())
# library(EasyMultiOmics)
library(openxlsx)
ps02 = EasyMultiOmics::ps.16s
ps01 = EasyMultiOmics::ps.its

#0 细菌/真菌数据合并----
ps03 =  merge16S_ITS(ps16s = ps01,
                     psITS = ps02,
                     NITS = 500,
                     N16s = 500,
                     scale = FALSE,
                     onlygroup =FALSE)

ps03

# 注意tax文件所产生的变化，OTU的名字加上了前缀，增加了filed列，定义细菌还是真菌数据
tax= ps03%>% vegan_tax() %>% as.data.frame()
head(tax)

# otu文件的变化
otu = ps03 %>% vegan_otu() %>% t() %>% as.data.frame()
head(otu)

#细菌真菌分组map文件要完全一致
map = ps03 %>% sample_data()
head(map)
map
# 定义后续参数-----

# 提取分组因子数量
gnum = phyloseq::sample_data(ps03)$Group %>% unique() %>% length()
gnum
#--设定排序顺序1：按照ps.16s对象中map文件顺序进行
axis_order =  phyloseq::sample_data(ps03)$Group %>%unique();axis_order

#-主题--颜色等
package.amp()
res = theme_my(ps03)
mytheme1 = res[[1]]
mytheme2 = res[[2]];
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]
col.g =  c("KO" = "#D55E00", "WT" = "#0072B2", "OE" = "#009E73")

# 创建细菌-真菌联合分析目录
bac_fun_path <- "./result/bacteria_fungi/"
dir.create(bac_fun_path, recursive = TRUE)

# alpha-diversity-----
#1 alpha.omics: 6中alpha多样性计算 ----
# 创建Alpha多样性分析目录
bac_fun_alpha_path <- file.path(bac_fun_path, "alpha_diversity")
dir.create(bac_fun_alpha_path, recursive = TRUE)

all.alpha = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )
tab = alpha.omics(ps = ps03,group = "Group")
head(tab)

data = cbind(data.frame(ID = 1:length(tab$Group),group = tab$Group),tab[all.alpha])
head(data)
data$ID = as.character(data$ID)

result = MuiKwWlx2(data = data,num =3:6)
result1 = FacetMuiPlotresultBox(data = data,num = 3:6,
                                result = result,
                                sig_show ="abc",ncol = 4 )
p1_1 = result1[[1]]
p1_1_final = p1_1+
  scale_fill_manual(values = col.g)+theme_classic()

res = FacetMuiPlotresultBar(data = data,num = c(3:(6)),result = result,sig_show ="abc",ncol = 4)
p1_2 = res[[1]]
p1_2_final = p1_2+
  scale_fill_manual(values = col.g)+theme_classic()

res = FacetMuiPlotReBoxBar(data = data,num = c(3:(6)),result = result,sig_show ="abc",ncol = 3)
p1_3 = res[[1]]
p1_3_final = p1_3+
  scale_fill_manual(values = col.g)+theme_classic()

#基于输出数据使用ggplot出图
p1_0 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) +
  geom_violin(alpha=1, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
  labs(x="", y="")+
  facet_wrap(.~name,scales="free_y",ncol  = 4) +
  # theme_classic()+
  geom_text(aes(x=group , y=y ,label=stat))
p1_0_final = p1_0+
  scale_fill_manual(values = col.g)+theme_classic()

# 保存Alpha多样性结果
ggsave(file.path(bac_fun_alpha_path, "alpha_diversity_box.png"), plot = p1_1_final, width = 12, height = 8, dpi = 300)
ggsave(file.path(bac_fun_alpha_path, "alpha_diversity_box.pdf"), plot = p1_1_final, width = 12, height = 8)
ggsave(file.path(bac_fun_alpha_path, "alpha_diversity_bar.png"), plot = p1_2_final, width = 12, height = 8, dpi = 300)
ggsave(file.path(bac_fun_alpha_path, "alpha_diversity_bar.pdf"), plot = p1_2_final, width = 12, height = 8)
ggsave(file.path(bac_fun_alpha_path, "alpha_diversity_boxbar.png"), plot = p1_3_final, width = 12, height = 8, dpi = 300)
ggsave(file.path(bac_fun_alpha_path, "alpha_diversity_boxbar.pdf"), plot = p1_3_final, width = 12, height = 8)
ggsave(file.path(bac_fun_alpha_path, "alpha_diversity_violin.png"), plot = p1_0_final, width = 12, height = 8, dpi = 300)
ggsave(file.path(bac_fun_alpha_path, "alpha_diversity_violin.pdf"), plot = p1_0_final, width = 12, height = 8)

# 创建Alpha多样性分析总表
bac_fun_alpha_wb <- createWorkbook()
addWorksheet(bac_fun_alpha_wb, "alpha_diversity_data")
writeData(bac_fun_alpha_wb, "alpha_diversity_data", tab, rowNames = TRUE)
addWorksheet(bac_fun_alpha_wb, "alpha_statistics")
writeData(bac_fun_alpha_wb, "alpha_statistics", result, rowNames = TRUE)
saveWorkbook(bac_fun_alpha_wb, file.path(bac_fun_alpha_path, "alpha_diversity.xlsx"), overwrite = TRUE)

# beta diversity-----
#2 ordinate.omics: 排序分析#----------
# 创建Beta多样性分析目录
bac_fun_beta_path <- file.path(bac_fun_path, "beta_diversity")
dir.create(bac_fun_beta_path, recursive = TRUE)

result = ordinate.omics(ps = ps03, group = "Group", dist = "bray",
                        method = "PCoA", Micromet = "anosim", pvalue.cutoff = 0.05,
                        pair = F)
p3_1 = result[[1]]
p3_1_final = p3_1+
  scale_fill_manual(values = col.g)+theme_classic()

#带标签图形出图
p3_2 = result[[3]]
p3_2_final = p3_2+
  scale_fill_manual(values = col.g)+theme_classic()

#---------排序-精修图
plotdata =result[[2]]
head(plotdata)
# 求均值
cent <- aggregate(cbind(x,y) ~Group, data = plotdata, FUN = mean)
cent
# 合并到样本坐标数据中
segs <- merge(plotdata, setNames(cent, c('Group','oNMDS1','oNMDS2')),
              by = 'Group', sort = FALSE)

library(ggsci)
p3_3 = p3_1 +geom_segment(data = segs,
                          mapping = aes(xend = oNMDS1, yend = oNMDS2,color = Group),show.legend=F) + # spiders
  geom_point(mapping = aes(x = x, y = y),data = cent, size = 5,pch = 24,color = "black",fill = "yellow")
p3_3_final = p3_3+  scale_fill_manual(values = col.g)+theme_classic()

# 保存Beta多样性结果
ggsave(file.path(bac_fun_beta_path, "PCoA_plot.png"), plot = p3_1_final, width = 8, height = 6, dpi = 300)
ggsave(file.path(bac_fun_beta_path, "PCoA_plot.pdf"), plot = p3_1_final, width = 8, height = 6)
ggsave(file.path(bac_fun_beta_path, "PCoA_with_labels.png"), plot = p3_2_final, width = 8, height = 6, dpi = 300)
ggsave(file.path(bac_fun_beta_path, "PCoA_with_labels.pdf"), plot = p3_2_final, width = 8, height = 6)
ggsave(file.path(bac_fun_beta_path, "PCoA_refined.png"), plot = p3_3_final, width = 8, height = 6, dpi = 300)
ggsave(file.path(bac_fun_beta_path, "PCoA_refined.pdf"), plot = p3_3_final, width = 8, height = 6)

# 创建Beta多样性分析总表
bac_fun_beta_wb <- createWorkbook()
addWorksheet(bac_fun_beta_wb, "PCoA_coordinates")
writeData(bac_fun_beta_wb, "PCoA_coordinates", plotdata, rowNames = TRUE)
addWorksheet(bac_fun_beta_wb, "centroids")
writeData(bac_fun_beta_wb, "centroids", cent, rowNames = TRUE)

#3 ordinateTest.omics:群落水平差异检测#-------
dat1 = ordinateTest.omics(ps = ps03, Micromet = "adonis", dist = "bray")
dat1

# 保存总体差异检测结果
addWorksheet(bac_fun_beta_wb, "overall_test")
writeData(bac_fun_beta_wb, "overall_test", dat1, rowNames = TRUE)

#4 pairordinateTest.omics:两两分组群落水平差异检测#-------
dat2 = pairordinateTest.omics(ps = ps03, Micromet = "MRPP", dist = "bray")
dat2

# 保存两两比较结果
addWorksheet(bac_fun_beta_wb, "pairwise_test")
writeData(bac_fun_beta_wb, "pairwise_test", dat2, rowNames = TRUE)

#5 mantal.omics ：群落差异检测普鲁士分析#------
dat = mantal.omics(ps01= ps01,ps02= ps02)

# 保存mantel检验结果
addWorksheet(bac_fun_beta_wb, "mantel_test")
writeData(bac_fun_beta_wb, "mantel_test", dat, rowNames = TRUE)

#6 cluster.omics:样品聚类#-----
res = cluster.omics (ps= ps03,hcluter_method = "complete",
                     dist = "bray",
                     cuttree = 3,
                     row_cluster = TRUE,
                     col_cluster =  TRUE)
p4 = res[[1]]
p4
p4_1 = res[[2]]
p4_1
p4_2 = res[[3]]
p4_2
dat = res[[4]]
dat

# 保存聚类分析结果
ggsave(file.path(bac_fun_beta_path, "cluster_plot1.png"), plot = p4, width = 10, height = 8, dpi = 300)
ggsave(file.path(bac_fun_beta_path, "cluster_plot1.pdf"), plot = p4, width = 10, height = 8)
ggsave(file.path(bac_fun_beta_path, "cluster_plot2.png"), plot = p4_1, width = 10, height = 8, dpi = 300)
ggsave(file.path(bac_fun_beta_path, "cluster_plot2.pdf"), plot = p4_1, width = 10, height = 8)
ggsave(file.path(bac_fun_beta_path, "cluster_plot3.png"), plot = p4_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(bac_fun_beta_path, "cluster_plot3.pdf"), plot = p4_2, width = 10, height = 8)
addWorksheet(bac_fun_beta_wb, "cluster_data")
writeData(bac_fun_beta_wb, "cluster_data", dat, rowNames = TRUE)
saveWorkbook(bac_fun_beta_wb, file.path(bac_fun_beta_path, "beta_diversity.xlsx"), overwrite = TRUE)

# constraint analysis----
# 创建约束分析目录
bac_fun_constraint_path <- file.path(bac_fun_path, "constraint_analysis")
dir.create(bac_fun_constraint_path, recursive = TRUE)

ps.tem = ps.16s %>% tax_glom_wt("Genus")
# 真菌筛选
res = loadingPCA.micro(ps = ps.its%>% tax_glom_wt("Genus"))
# p = res[[1]]
# p
dat = res[[2]]
dat$id = row.names(dat)
id = dat %>% arrange(desc(PCone)) %>% head(15) %>% .$id
id

ftab = ps.its%>% tax_glom_wt("Genus") %>% scale_micro() %>%
  phyloseq::subset_taxa(
    row.names(phyloseq::tax_table(ps.its%>% tax_glom_wt("Genus"))) %in% c(id)
  ) %>% vegan_otu() %>%
  as.data.frame()

head(ftab)

#7 RDA_CCA:细菌真菌 共排序-----
result = RDA_CCA(ps = ps.tem,
                 env = ftab,
                 # path = RDApath,
                 chose.env = FALSE
)

p1 = result[[1]]
p1_final = p1+
  scale_fill_manual(values = col.g)+theme_classic()

# 提取作图数据
dataplot = result[[2]]
dataplot
# 提取带有标记的图片
p2 = result[[3]]
p2_final = p2+
  scale_fill_manual(values = col.g)+theme_classic()
aov = result[[4]]
aov

# 保存RDA/CCA结果
ggsave(file.path(bac_fun_constraint_path, "RDA_CCA_plot1.png"), plot = p1_final, width = 8, height = 6, dpi = 300)
ggsave(file.path(bac_fun_constraint_path, "RDA_CCA_plot1.pdf"), plot = p1_final, width = 8, height = 6)
ggsave(file.path(bac_fun_constraint_path, "RDA_CCA_plot2.png"), plot = p2_final, width = 8, height = 6, dpi = 300)
ggsave(file.path(bac_fun_constraint_path, "RDA_CCA_plot2.pdf"), plot = p2_final, width = 8, height = 6)

# 创建约束分析总表
bac_fun_constraint_wb <- createWorkbook()
addWorksheet(bac_fun_constraint_wb, "RDA_CCA_data")
writeData(bac_fun_constraint_wb, "RDA_CCA_data", dataplot, rowNames = TRUE)
addWorksheet(bac_fun_constraint_wb, "RDA_CCA_anova")
writeData(bac_fun_constraint_wb, "RDA_CCA_anova", aov, rowNames = TRUE)
addWorksheet(bac_fun_constraint_wb, "fungi_selected")
writeData(bac_fun_constraint_wb, "fungi_selected", ftab, rowNames = TRUE)

#8 RDA_CCA_explain_percent: 对微生物群落的解释比例-----
result = RDA_CCA_explain_percent(ps = ps.tem, env.dat = ftab)

sam = result[[1]]
sam

all = result[[2]]

p= result[[3]]
p

# 保存解释比例结果
ggsave(file.path(bac_fun_constraint_path, "explain_percent_plot.png"), plot = p, width = 8, height = 6, dpi = 300)
ggsave(file.path(bac_fun_constraint_path, "explain_percent_plot.pdf"), plot = p, width = 8, height = 6)
addWorksheet(bac_fun_constraint_wb, "explain_sample")
writeData(bac_fun_constraint_wb, "explain_sample", sam, rowNames = TRUE)
addWorksheet(bac_fun_constraint_wb, "explain_all")
writeData(bac_fun_constraint_wb, "explain_all", all, rowNames = TRUE)

#9 rdacca.hp.micro:  层次分割对细菌群落影响的真菌---------
library(rdacca.hp)
library(vegan)
library(ggClusterNet)
library(tidyverse)
res = rdacca.hp.micro(OTU = ps.tem %>% filter_OTU_ps(500) %>%vegan_otu() %>% as.data.frame(), env = ftab[,1:5], cca = FALSE)
p3_rda = res[[1]]
p3_rda_final = p3_rda+xlab("")
dat1 = res[[2]]
dat1
p3_db_rda  = res[[3]]
p3_db_rda_final = p3_db_rda+xlab("")
dat2 = res[[4]]
dat2

# 保存层次分割结果
ggsave(file.path(bac_fun_constraint_path, "hierarchical_partition1.png"), plot = p3_rda_final, width = 8, height = 6, dpi = 300)
ggsave(file.path(bac_fun_constraint_path, "hierarchical_partition1.pdf"), plot = p3_rda_final, width = 8, height = 6)
ggsave(file.path(bac_fun_constraint_path, "hierarchical_partition2.png"), plot = p3_db_rda_final, width = 8, height = 6, dpi = 300)
ggsave(file.path(bac_fun_constraint_path, "hierarchical_partition2.pdf"), plot = p3_db_rda_final, width = 8, height = 6)
addWorksheet(bac_fun_constraint_wb, "hierarchical_data1")
writeData(bac_fun_constraint_wb, "hierarchical_data1", dat1, rowNames = TRUE)
addWorksheet(bac_fun_constraint_wb, "hierarchical_data2")
writeData(bac_fun_constraint_wb, "hierarchical_data2", dat2, rowNames = TRUE)

#10 rdacca.hp.micro.p: 层次分割对微生物影响的代谢物#---------
res <-rdacca.hp.micro.p(
  OTU = ps.tem %>% filter_OTU_ps(200) %>%vegan_otu() %>% as.data.frame(),
  env = ftab[,1:5],
  cca = FALSE,
  dbRDA = FALSE
)
p3_rda = res[[1]]
p3_rda_final2 = p3_rda+xlab("")
dat1 = res[[2]]
dat1
p3_db_rda  = res[[3]]
p3_db_rda_final2 = p3_db_rda+xlab("")

# 保存带p值的层次分割结果
ggsave(file.path(bac_fun_constraint_path, "hierarchical_partition_p1.png"), plot = p3_rda_final2, width = 8, height = 6, dpi = 300)
ggsave(file.path(bac_fun_constraint_path, "hierarchical_partition_p1.pdf"), plot = p3_rda_final2, width = 8, height = 6)
ggsave(file.path(bac_fun_constraint_path, "hierarchical_partition_p2.png"), plot = p3_db_rda_final2, width = 8, height = 6, dpi = 300)
ggsave(file.path(bac_fun_constraint_path, "hierarchical_partition_p2.pdf"), plot = p3_db_rda_final2, width = 8, height = 6)
addWorksheet(bac_fun_constraint_wb, "hierarchical_p_data")
writeData(bac_fun_constraint_wb, "hierarchical_p_data", dat1, rowNames = TRUE)
saveWorkbook(bac_fun_constraint_wb, file.path(bac_fun_constraint_path, "constraint_analysis.xlsx"), overwrite = TRUE)

# 网络分析---------
# 创建网络分析目录
bac_fun_network_path <- file.path(bac_fun_path, "network_analysis")
dir.create(bac_fun_network_path, recursive = TRUE)

#11 corBionetwork.st-----
pst = ps03 %>% tax_glom_wt("Genus")
res = corBionetwork.st(
  ps.st= ps03,# phyloseq对象
  N= 500,
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",
  method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,
  layout_net = "model_maptree2",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 5,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,
  maxsize = 14)

dat = res[[2]]

cor = dat$cortab
p = res[[1]]
p

# 保存网络图
ggsave(file.path(bac_fun_network_path, "network_plot.png"), plot = p, width = 12, height = 10, dpi = 300)
ggsave(file.path(bac_fun_network_path, "network_plot.pdf"), plot = p, width = 12, height = 10)

# 创建网络分析总表
bac_fun_network_wb <- createWorkbook()

#12 net_properties.4-------
id = names(cor)
id

for (i in 1:length(id)) {
  igraph= cor[[id[i]]] %>% make_igraph()
  dat = net_properties.4(igraph,n.hub = FALSE)
  head(dat,n = 16)
  colnames(dat) = id[i]

  if (i == 1) {
    dat2 = dat
  } else{
    dat2 = cbind(dat2,dat)
  }
}
head(dat2)

# 保存网络属性结果
addWorksheet(bac_fun_network_wb, "network_properties")
writeData(bac_fun_network_wb, "network_properties", dat2, rowNames = TRUE)

#13 netproperties.sample-----
id = pst %>% sample_data() %>% .$Group %>% unique()
id
id0 = names(cor)
for (i in 1:length(id)) {
  pst1 = pst %>% subset_samples.wt("Group",id[i]) %>% remove.zero()
  dat.f = netproperties.sample(pst = pst1,cor = cor[[id0[i]]])
  # head(dat.f)
  if (i == 1) {
    dat.f2 = dat.f
  } else{
    dat.f2 = rbind(dat.f2,dat.f)
  }
}

# 保存样本网络属性结果
addWorksheet(bac_fun_network_wb, "sample_network_properties")
writeData(bac_fun_network_wb, "sample_network_properties", dat.f2, rowNames = TRUE)

#31 module.compare.net.pip -----

dat = module.compare.net.pip(
  ps = NULL,
  corg = cor,
  degree = TRUE,
  zipi = FALSE,
  r.threshold= 0.8,
  p.threshold=0.05,
  method = "spearman",
  padj = F,
  n = 3)

res = dat[[1]]

# 保存网络模块比较结果
addWorksheet(bac_fun_network_wb, "module_comparison")
writeData(bac_fun_network_wb, "module_comparison", res, rowNames = TRUE)

#14 node_properties-------

for (i in 1:length(id0)) {

  igraph= cor[[id0[i]]] %>% make_igraph()
  nodepro = node_properties(igraph) %>% as.data.frame()
  nodepro$Group = id[i]
  head(nodepro)
  colnames(nodepro) = paste0(colnames(nodepro),".",id[i])
  nodepro = nodepro %>%
    as.data.frame() %>%
    rownames_to_column("ASV.name")

  # head(dat.f)
  if (i == 1) {
    nodepro2 = nodepro
  } else{

    nodepro2 = nodepro2 %>% full_join(nodepro,by = "ASV.name")

  }
}
head(nodepro2)

# 保存节点属性结果
addWorksheet(bac_fun_network_wb, "node_properties")
writeData(bac_fun_network_wb, "node_properties", nodepro2, rowNames = TRUE)
saveWorkbook(bac_fun_network_wb, file.path(bac_fun_network_path, "network_analysis.xlsx"), overwrite = TRUE)

# 机器学习----
# 创建机器学习目录
bac_fun_ml_path <- file.path(bac_fun_path, "machine_learning")
dir.create(bac_fun_ml_path, recursive = TRUE)

library(randomForest)
library(caret)
library(ROCR)
library(e1071)
library(rpart)
library(mia)
library(TreeSummarizedExperiment)
library("class")
library(kknn)
library(ipred)
library(glmnet)
library(nnet)
library(MLeval)

#rfcv.omics :交叉验证结果-------
ps = ps03 %>% tax_glom_wt("Genus")

result =rfcv.omics(ps = ps ,group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
prfcv_final = prfcv+theme_classic()

rfcvtable = result[[3]]
rfcvtable

# 保存交叉验证结果
ggsave(file.path(bac_fun_ml_path, "rfcv_plot.png"), plot = prfcv_final, width = 8, height = 6, dpi = 300)
ggsave(file.path(bac_fun_ml_path, "rfcv_plot.pdf"), plot = prfcv_final, width = 8, height = 6)

# 创建机器学习分析总表
bac_fun_ml_wb <- createWorkbook()
addWorksheet(bac_fun_ml_wb, "rfcv_results")
writeData(bac_fun_ml_wb, "rfcv_results", rfcvtable, rowNames = TRUE)

#8 Roc.omics:ROC 曲线绘制----
id = sample_data(ps03)$Group %>% unique()
aaa = combn(id,2)
i= 1
group = c(aaa[1,i],aaa[2,i])
b= data.frame(group)
sample_data(ps03)

pst = ps03 %>% subset_samples.wt("Group",group)

res = Roc.omics( ps = pst %>% filter_OTU_ps(100),repnum = 5)
p33.1 =  res[[1]]
p33.1_final = p33.1+theme_classic()
p33.2 =  res[[2]]
p33.2
dat =  res[[3]]
dat

# 保存ROC曲线结果
ggsave(file.path(bac_fun_ml_path, "ROC_plot1.png"), plot = p33.1_final, width = 8, height = 6, dpi = 300)
ggsave(file.path(bac_fun_ml_path, "ROC_plot1.pdf"), plot = p33.1_final, width = 8, height = 6)
ggsave(file.path(bac_fun_ml_path, "ROC_plot2.png"), plot = p33.2, width = 8, height = 6, dpi = 300)
ggsave(file.path(bac_fun_ml_path, "ROC_plot2.pdf"), plot = p33.2, width = 8, height = 6)
addWorksheet(bac_fun_ml_wb, "ROC_data")
writeData(bac_fun_ml_wb, "ROC_data", dat, rowNames = TRUE)

#9 loadingPCA.omics:载荷矩阵筛选特征功能------
res = loadingPCA.omics(ps = ps03,Top = 20)
p34.1 = res[[1]]
p34.1_final = p34.1+theme_classic()+ ylab("")
dat = res[[2]]
dat

# 保存PCA载荷结果
ggsave(file.path(bac_fun_ml_path, "PCA_loading_plot.png"), plot = p34.1_final, width = 10, height = 8, dpi = 300)
ggsave(file.path(bac_fun_ml_path, "PCA_loading_plot.pdf"), plot = p34.1_final, width = 10, height = 8)
addWorksheet(bac_fun_ml_wb, "PCA_loading")
writeData(bac_fun_ml_wb, "PCA_loading", dat, rowNames = TRUE)

#10 randomforest.omics: 随机森林筛选特征功能----
res = randomforest.omics( ps = ps03 %>% filter_OTU_ps(100),
                          group  = "Group",
                          optimal = 50)
p42.1 = res[[1]]
p42.1_final = p42.1+theme_classic()+ylab("")
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
p42.4 =res[[4]]
p42.4

# 保存随机森林结果
ggsave(file.path(bac_fun_ml_path, "randomforest_plot1.png"), plot = p42.1_final, width = 10, height = 8, dpi = 300)
ggsave(file.path(bac_fun_ml_path, "randomforest_plot1.pdf"), plot = p42.1_final, width = 10, height = 8)
ggsave(file.path(bac_fun_ml_path, "randomforest_plot2.png"), plot = p42.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(bac_fun_ml_path, "randomforest_plot2.pdf"), plot = p42.2, width = 10, height = 8)
ggsave(file.path(bac_fun_ml_path, "randomforest_plot4.png"), plot = p42.4, width = 10, height = 8, dpi = 300)
ggsave(file.path(bac_fun_ml_path, "randomforest_plot4.pdf"), plot = p42.4, width = 10, height = 8)
addWorksheet(bac_fun_ml_wb, "randomforest_features")
writeData(bac_fun_ml_wb, "randomforest_features", dat, rowNames = TRUE)

#11 svm: 筛选特征细菌与真菌------
res <- svm.omics(ps = pst %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存SVM结果
addWorksheet(bac_fun_ml_wb, "SVM_AUC")
writeData(bac_fun_ml_wb, "SVM_AUC", AUC, rowNames = TRUE)
addWorksheet(bac_fun_ml_wb, "SVM_importance")
writeData(bac_fun_ml_wb, "SVM_importance", importance, rowNames = TRUE)

#12 GLM: 筛选特征细菌与真菌---------
res <- glm.omics(ps = pst %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存GLM结果
addWorksheet(bac_fun_ml_wb, "GLM_AUC")
writeData(bac_fun_ml_wb, "GLM_AUC", AUC, rowNames = TRUE)
addWorksheet(bac_fun_ml_wb, "GLM_importance")
writeData(bac_fun_ml_wb, "GLM_importance", importance, rowNames = TRUE)

#13 xgboost: 筛选特征细菌与真菌-----
res =xgboost.omics(ps = ps03, top = 5,k =5  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存XGBoost结果
addWorksheet(bac_fun_ml_wb, "XGBoost_accuracy")
writeData(bac_fun_ml_wb, "XGBoost_accuracy", accuracy, rowNames = TRUE)
addWorksheet(bac_fun_ml_wb, "XGBoost_importance")
writeData(bac_fun_ml_wb, "XGBoost_importance", importance, rowNames = TRUE)

#14 kNN: 筛选特征细菌与真菌--------
res = knn.omics(ps =ps03, top = 20,k =5 )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存kNN结果
addWorksheet(bac_fun_ml_wb, "kNN_accuracy")
writeData(bac_fun_ml_wb, "kNN_accuracy", accuracy, rowNames = TRUE)
addWorksheet(bac_fun_ml_wb, "kNN_importance")
writeData(bac_fun_ml_wb, "kNN_importance", importance, rowNames = TRUE)

# #15 decisiontree: 筛选特征细菌与真菌------
res =decisiontree.omics(ps=ps03, top =500,  k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存决策树结果
addWorksheet(bac_fun_ml_wb, "DecisionTree_accuracy")
writeData(bac_fun_ml_wb, "DecisionTree_accuracy", accuracy, rowNames = TRUE)
addWorksheet(bac_fun_ml_wb, "DecisionTree_importance")
writeData(bac_fun_ml_wb, "DecisionTree_importance", importance, rowNames = TRUE)

#16 bagging: 筛选特征细菌与真菌-----
res =bagging.omics(ps =  ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Bagging结果
addWorksheet(bac_fun_ml_wb, "Bagging_accuracy")
writeData(bac_fun_ml_wb, "Bagging_accuracy", accuracy, rowNames = TRUE)
addWorksheet(bac_fun_ml_wb, "Bagging_importance")
writeData(bac_fun_ml_wb, "Bagging_importance", importance, rowNames = TRUE)

#17 lasso: 筛选特征细菌与真菌-----
res = lasso.omics (ps =  ps03, top = 20, seed = 10, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Lasso结果
addWorksheet(bac_fun_ml_wb, "Lasso_accuracy")
writeData(bac_fun_ml_wb, "Lasso_accuracy", accuracy, rowNames = TRUE)
addWorksheet(bac_fun_ml_wb, "Lasso_importance")
writeData(bac_fun_ml_wb, "Lasso_importance", importance, rowNames = TRUE)

#18 naivebayes: 筛选特征细菌与真菌------
res = naivebayes.omics(ps=ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存朴素贝叶斯结果
addWorksheet(bac_fun_ml_wb, "NaiveBayes_accuracy")
writeData(bac_fun_ml_wb, "NaiveBayes_accuracy", accuracy, rowNames = TRUE)
addWorksheet(bac_fun_ml_wb, "NaiveBayes_importance")
writeData(bac_fun_ml_wb, "NaiveBayes_importance", importance, rowNames = TRUE)

# #19 nnet神经网络: 筛选特征细菌与真菌------
res =nnet.omics(ps=ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存神经网络结果
addWorksheet(bac_fun_ml_wb, "NeuralNet_accuracy")
writeData(bac_fun_ml_wb, "NeuralNet_accuracy", accuracy, rowNames = TRUE)
addWorksheet(bac_fun_ml_wb, "NeuralNet_importance")
writeData(bac_fun_ml_wb, "NeuralNet_importance", importance, rowNames = TRUE)

# 保存机器学习分析总表
saveWorkbook(bac_fun_ml_wb, file.path(bac_fun_ml_path, "machine_learning.xlsx"), overwrite = TRUE)
