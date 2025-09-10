#清空内存#######
rm(list=ls())
# R包加载
library(sva)          # ComBat批次校正
library(RSpectra)     # 稀疏PCA
library(dplyr)        # 数据操作
library(compositions) # CLR
library(ggClusterNet)
library(phyloseq)
library(devtools)
load_all()
library(EasyMultiOmics)
library(caret)
library("DESeq2")
library(limma)
library(randomForest)
library(DESeq2)
library(igraph)
library(openxlsx)

# 创建主要目录结构 ----
main_dir <- "./result/across_omics"
if (!dir.exists(main_dir)) dir.create(main_dir)


# 挑选基因组学（以微生物组为例）中的特征#------
ps = ps.16s %>% subset_samples.wt("Group",c("WT","KO"))
# ps = ps %>% subset_taxa.wt("Species",c("unclassified"),TRUE)
ps = ps %>% tax_glom_wt(6)
ps

res = count_selection(ps,
                      group_var = "Group",
                      prevalence_threshold = 0.1,
                      detection_threshold = 0.001,
                      diff_method = "DESeq2",
                      ml_method = "rf",
                      cor_method = "spearman",
                      cor_threshold = 0.6,
                      p_threshold = 0.05,
                      weights = list(abundance = 0.3,
                                     importance = 0.3,
                                     differential = 0.2,
                                     network = 0.2))

dat = res$final_ranking %>% head(100)

# 保存基因组特征选择结果
genomics_wb <- createWorkbook()
addWorksheet(genomics_wb, "genomics_features")
addWorksheet(genomics_wb, "genomics_selection_stats")
writeData(genomics_wb, "genomics_features", dat, rowNames = TRUE)
if(!is.null(res$stats)) {
  writeData(genomics_wb, "genomics_selection_stats", res$stats, rowNames = TRUE)
}
saveWorkbook(genomics_wb, file.path(main_dir, "genomics_feature_selection_results.xlsx"), overwrite = TRUE)


id = dat$OTU[1]
id

##寻找与特征微生物潜在相关的代谢物------
ps2 = ps.ms %>% subset_samples.wt("Group",c("WT","OE"))
tab = ps2 %>% vegan_otu() %>%
  as.data.frame()
ps = ps.16s %>% subset_samples.wt("Group",c("WT","OE")) %>% tax_glom_wt(6)

# 寻找与特征微生物相关联的代谢物 Find Associated Metabolites with Target Microbe

?count_to_mass
results <- count_to_mass(
  phyloseq_obj = ps,
  metabolome_mat = tab,
  target_microbe_name = "Actinocorallia" ,
  top_n = 30
)

names(results)
results$top_metabolites
tax = ps2 %>% tax_table() %>% as.data.frame()
head(tax)

dat = results$top_metabolites %>% left_join(tax,by = c("metabolite"="metab_id"))

head(dat)
dat$Metabolite
# 查看结果
print(head(results$top_metabolites))

# 保存微生物-代谢物关联结果
microbe_metabolite_wb <- createWorkbook()
addWorksheet(microbe_metabolite_wb, "microbe_metabolite_association")
addWorksheet(microbe_metabolite_wb, "association_details")
writeData(microbe_metabolite_wb, "microbe_metabolite_association", dat, rowNames = TRUE)
if(!is.null(results$correlation_matrix)) {
  writeData(microbe_metabolite_wb, "association_details", results$correlation_matrix, rowNames = TRUE)
}
saveWorkbook(microbe_metabolite_wb, file.path(main_dir, "microbe_metabolite_association_results.xlsx"), overwrite = TRUE)

#  multi_omics_alignment 多组学矫正#-----
# 对输入数据进行校正后，寻找组学特征以及潜在关联
otu = ps.16s %>%
  tax_glom_wt(6) %>%
  vegan_otu() %>%
  as.data.frame()

tab = ps.ms %>% vegan_otu() %>%
  as.data.frame()

omics_list = list(micro = otu,ms = tab)

# group = sample_data(ps.16s)$Group
# names(group) = row.names(sample_data(ps.16s))
# microbiome_data = otu  # 样本×微生物
# metabolome_data = tab  # 样本×代谢物
# n_factors = 10                   # 潜在因子数
# method = "auto"

alignment_results <- multi_omics_alignment(
  microbiome_data = otu,  # 样本×微生物
  metabolome_data = tab,  # 样本×代谢物
  n_factors = 10,                       # 潜在因子数
  method = "fast"                       # 自动选择快速或完整方法
)

# 查看对齐质量
print(alignment_results$quality_details)

alignment_results$aligned_microbiome
alignment_results$aligned_metabolome

# 保存多组学对齐结果
alignment_wb <- createWorkbook()
addWorksheet(alignment_wb, "alignment_quality")
addWorksheet(alignment_wb, "aligned_microbiome")
addWorksheet(alignment_wb, "aligned_metabolome")
writeData(alignment_wb, "alignment_quality", alignment_results$quality_details, rowNames = TRUE)
writeData(alignment_wb, "aligned_microbiome", alignment_results$aligned_microbiome, rowNames = TRUE)
writeData(alignment_wb, "aligned_metabolome", alignment_results$aligned_metabolome, rowNames = TRUE)
saveWorkbook(alignment_wb, file.path(main_dir, "multi_omics_alignment_results.xlsx"), overwrite = TRUE)


ps1 = ps.16s %>% tax_glom_wt(6)
otu_table(ps1) = otu_table(as.matrix(alignment_results$aligned_microbiome) %>% t(),taxa_are_rows = TRUE)

ps2 = ps.ms
otu_table(ps2) = otu_table(as.matrix(alignment_results$aligned_metabolome) %>% t(),taxa_are_rows = TRUE)

ps2 = ps2 %>% subset_samples.wt("Group",c("WT","KO"))
tab = ps2 %>% vegan_otu() %>%
  as.data.frame()


ps = ps1 %>% subset_samples.wt("Group",c("WT","KO"))

# phyloseq_obj = ps
# metabolome_mat = tab
# target_microbe_name = "Arthrobacter"
# top_n = 30


## find_associated_metabolites.2 计算关联#-----
# 运行完整分析
results <- count_to_mass(
  phyloseq_obj = ps,
  metabolome_mat = tab,
  target_microbe_name = "Arthrobacter",
  top_n = 50,
  use_normalization = FALSE
)

names(results)
results$top_metabolites

tax = ps2 %>% tax_table() %>% as.data.frame()
head(tax)
dat = results$top_metabolites %>% left_join(tax,by = c("metabolite"="metab_id"))
head(dat)
dat %>% tail()
dat$Metabolite

# 保存关联代谢物分析结果
associated_metabolites_wb <- createWorkbook()
addWorksheet(associated_metabolites_wb, "associated_metabolites")
addWorksheet(associated_metabolites_wb, "correlation_data")
writeData(associated_metabolites_wb, "associated_metabolites", dat, rowNames = TRUE)
if(!is.null(results$correlation_data)) {
  writeData(associated_metabolites_wb, "correlation_data", results$correlation_data, rowNames = TRUE)
}
saveWorkbook(associated_metabolites_wb, file.path(main_dir, "associated_metabolites_analysis_results.xlsx"), overwrite = TRUE)

## smart_metabolite_selection 指定关联代谢物数量#------
selection_report <- feature_selection(
  results = results,
  sample_size = 20,           # 你的样本量
  validation_available = F, # 是否有验证队列
  resource_level = "low"    # 资源水平: "low", "medium", "high"
)

ggsave(file.path(main_dir, "method_comparison.pdf"), plot = selection_report$plots[[1]], width = 12, height = 8, dpi = 300)
ggsave(file.path(main_dir, "method_comparison.png"), plot = selection_report$plots[[1]], width = 12, height = 8, dpi = 300)
ggsave(file.path(main_dir, "score_distribution.pdf"), plot = selection_report$plots[[2]], width = 12, height = 8, dpi = 300)
ggsave(file.path(main_dir, "score_distribution.png"), plot = selection_report$plots[[2]], width = 12, height = 8, dpi = 300)

# 获取推荐的代谢物列表
selected_metabolites <- selection_report$metabolites$standard

tax = ps2 %>% tax_table() %>% as.data.frame()
head(tax)

dat = selected_metabolites %>% left_join(tax,by = c("metabolite"="metab_id"))
dat$Metabolite
head(dat)

# 保存智能代谢物选择结果
smart_selection_wb <- createWorkbook()
addWorksheet(smart_selection_wb, "smart_metabolite_selection")
addWorksheet(smart_selection_wb, "selection_report")
writeData(smart_selection_wb, "smart_metabolite_selection", dat, rowNames = TRUE)
if(!is.null(selection_report$report)) {
  writeData(smart_selection_wb, "selection_report", selection_report$report, rowNames = TRUE)
}
saveWorkbook(smart_selection_wb, file.path(main_dir, "smart_metabolite_selection_results.xlsx"), overwrite = TRUE)


#  代谢组/蛋白组等学特征挑选#-----
tab = ps.ms %>%subset_samples.wt("Group",c("WT","KO")) %>% vegan_otu( ) %>%
  as.data.frame()
?mass_selection
map = sample_data(ps.ms %>%subset_samples.wt("Group",c("WT","KO")) )
map$Group = as.factor(map$Group)
map$ID = NULL
## 2) 运行特征挑选主流程
res <- mass_selection(
  ms_mat         = tab,
  sample_meta    = map,
  group_var      = "Group",
  batch_var      = NULL,
  qc_flag_var    = NULL,
  prevalence_threshold = 0.7,      # ≥70% 样本非缺失才保留
  impute_method  = "halfmin",      # 缺失填补：半最小值
  normalize_method = "median",     # 归一化：样本中位数
  combat         = F,           # 批次校正：ComBat
  qc_drift_correction = TRUE,      # 用 QC 做 LOESS 漂移校正
  qc_rsd_cut     = 0.30,           # QC RSD ≤ 30% 保留
  diff_method    = "limma",        # 差异分析：limma（两组更稳）
  ml_method      = "rf",           # 机器学习：随机森林
  cor_method     = "spearman",     # 网络相关：斯皮尔曼
  cor_threshold  = 0.6,
  p_threshold    = 0.05,
  weights = list(
    abundance   = 0.3,   # 偏好高丰度与高流行率
    importance  = 0.3,   # ML 重要性
    differential= 0.2,   # 差异显著性与效应量
    network     = 0.2    # 网络中心性
  ),
  nfolds = 5,
  ntree  = 300           # RF 森林数（示例里稍微降点速度更快）
)

## 3) 查看结果
names(res)
head(res$final_ranking, 15)     # 综合排名前15

dat = res$final_ranking
head(dat)

tax = ps.ms %>% tax_table() %>% as.data.frame()
head(tax)

dat =  res$final_ranking %>% left_join(tax,by = c("Metabolite"="metab_id"))

head(dat)
dat$Metabolite.y %>% head(50)

# 保存代谢组特征选择结果
metabolomics_wb <- createWorkbook()
addWorksheet(metabolomics_wb, "metabolomics_features")
addWorksheet(metabolomics_wb, "selection_statistics")
writeData(metabolomics_wb, "metabolomics_features", dat, rowNames = TRUE)
if(!is.null(res$statistics)) {
  writeData(metabolomics_wb, "selection_statistics", res$statistics, rowNames = TRUE)
}
saveWorkbook(metabolomics_wb, file.path(main_dir, "metabolomics_feature_selection_results.xlsx"), overwrite = TRUE)

##----寻找与代谢物相关的微生物----
ps2 = ps.ms %>% subset_samples.wt("Group",c("WT","OE"))
tab = ps2 %>% vegan_otu() %>%
  as.data.frame()

ps = ps.16s %>% subset_samples.wt("Group",c("WT","OE")) %>% tax_glom_wt(6)

res_sim <- mass_to_count(
  phyloseq_obj = ps,
  metabolome_mat = tab,
  target_metabolite_name = "metab_7974",  # 我们植入信号的代谢物
  top_n = 20,
  cor_method = "spearman",
  prevalence_threshold = 0.10,
  detection_threshold  = 1e-4,
  use_clr = TRUE,
  group_split = "median"
)

res_sim$top_microbes$microbe

# 保存代谢物-微生物关联结果
metabolite_microbe_wb <- createWorkbook()
addWorksheet(metabolite_microbe_wb, "metabolite_microbe_association")
addWorksheet(metabolite_microbe_wb, "top_microbes_details")
writeData(metabolite_microbe_wb, "metabolite_microbe_association", res_sim$top_microbes, rowNames = TRUE)
if(!is.null(res_sim$correlation_matrix)) {
  writeData(metabolite_microbe_wb, "top_microbes_details", res_sim$correlation_matrix, rowNames = TRUE)
}
saveWorkbook(metabolite_microbe_wb, file.path(main_dir, "metabolite_microbe_association_results.xlsx"), overwrite = TRUE)




