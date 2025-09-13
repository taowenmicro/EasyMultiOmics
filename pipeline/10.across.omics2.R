# 清空内存
rm(list = ls())

# R包加载

library(RSpectra)     # 稀疏PCA
library(compositions) # CLR
library(ggClusterNet)
library(phyloseq)
library(EasyMultiOmics)
library(caret)
library("DESeq2")
library(limma)
library(randomForest)
library(DESeq2)
library(igraph)

library(patchwork)
library(scales)
library(tidygraph)
library(ggraph)
library(gridExtra)
library(openxlsx)

# 创建主要目录结构
main_dir <- "../result/across_omics"
if (!dir.exists(main_dir)) dir.create(main_dir, recursive = TRUE)

# 挑选基因组学中的特征 ------
cat("开始基因组学特征选择...\n")
ps = ps.16s %>% subset_samples.wt("Group",c("WT","KO"))
# ps = ps %>% subset_taxa.wt("Species",c("unclassified"),TRUE)
ps = ps %>% tax_glom_wt(6)
ps
detach(package:MicrobiotaProcess)
res = microbiome_feature_selection(ps,
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
addWorksheet(genomics_wb, "diff_results")
addWorksheet(genomics_wb, "importance_results")
addWorksheet(genomics_wb, "network_results")
addWorksheet(genomics_wb, "abundance_results")

writeData(genomics_wb, "genomics_features", dat, rowNames = TRUE)
writeData(genomics_wb, "genomics_selection_stats", res$stats, rowNames = TRUE)
writeData(genomics_wb, "diff_results", res$diff_results, rowNames = TRUE)
writeData(genomics_wb, "importance_results", res$importance_results, rowNames = TRUE)
writeData(genomics_wb, "network_results", res$network_results, rowNames = TRUE)
writeData(genomics_wb, "abundance_results", res$abundance_results, rowNames = TRUE)
saveWorkbook(genomics_wb, file.path(main_dir, "genomics_feature_selection_results.xlsx"), overwrite = TRUE)

id = dat$OTU[1]

# 绘制放射状图
p1 <- plot_radial_top_otus(dat, top_n=50, inner_ratio = 0.45)
print(p1)
ggsave(file.path(main_dir, "radial_top_otus.pdf"), p1, width = 10, height = 10)
ggsave(file.path(main_dir, "radial_top_otus.png"), p1, width = 10, height = 10, dpi = 300)

p2 <- plot_physeq_features(ps, dat, top_n=50)
print(p2)
ggsave(file.path(main_dir, "physeq_features.pdf"), p2, width = 12, height = 8)
ggsave(file.path(main_dir, "physeq_features.png"), p2, width = 12, height = 8, dpi = 300)


# 跨组学特征寻找 -------

ps2 = ps.ms %>% subset_samples.wt("Group",c("WT","KO"))
tab = ps2 %>% vegan_otu() %>%
  as.data.frame()

ps = ps.16s %>% subset_samples.wt("Group",c("WT","KO")) %>% tax_glom_wt(6)

# 运行完整分析
results <- find_associated_metabolites(
  phyloseq_obj = ps,
  metabolome_mat = tab,
  target_microbe_name = "Pseudonocardia",
  top_n = 30
)

names(results)
results$top_metabolites

tax = ps2 %>% tax_table() %>% as.data.frame()
head(tax)

dat = results$top_metabolites %>% left_join(tax, by = c("metabolite"="metab_id"))
head(dat)

# 保存微生物-代谢物关联结果
microbe_metabolite_wb <- createWorkbook()
addWorksheet(microbe_metabolite_wb, "microbe_metabolite_association")
addWorksheet(microbe_metabolite_wb, "association_details")
addWorksheet(microbe_metabolite_wb, "correlation_results")
writeData(microbe_metabolite_wb, "microbe_metabolite_association", dat, rowNames = TRUE)
writeData(microbe_metabolite_wb, "association_details", results$correlation_matrix, rowNames = TRUE)
writeData(microbe_metabolite_wb, "correlation_results", results$correlation_results, rowNames = TRUE)
saveWorkbook(microbe_metabolite_wb, file.path(main_dir, "microbe_metabolite_association_results.xlsx"), overwrite = TRUE)

p3 <- plot_microbe_circular_network2(dat,"Microbe_x",thr=0.5)
p3
ggsave(file.path(main_dir, "cross_omics_network.pdf"), p3, width = 10, height = 10)
ggsave(file.path(main_dir, "cross_omics_network.png"), p3, width = 10, height = 10, dpi = 300)

# multi_omics_alignment 多组学矫正 -----

otu = ps.16s %>%
  tax_glom_wt(6) %>%
  vegan_otu() %>%
  as.data.frame()

tab_ms = ps.ms %>% vegan_otu() %>%
  as.data.frame()

omics_list = list(micro = otu, ms = tab_ms)

# 基本使用
alignment_results <- multi_omics_alignment(
  microbiome_data = otu,
  metabolome_data = tab_ms,
  n_factors = 10,
  method = "fast"
)

# 查看对齐质量
print(alignment_results$quality_details)

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
otu_table(ps1) = otu_table(as.matrix(alignment_results$aligned_microbiome) %>% t(), taxa_are_rows = TRUE)

ps2_aligned = ps.ms
otu_table(ps2_aligned) = otu_table(as.matrix(alignment_results$aligned_metabolome) %>% t(), taxa_are_rows = TRUE)

ps2_aligned = ps2_aligned %>% subset_samples.wt("Group",c("WT","KO"))
tab_aligned = ps2_aligned %>% vegan_otu() %>%
  as.data.frame()

ps_aligned = ps1 %>% subset_samples.wt("Group",c("WT","KO"))

# find_associated_metabolites.2 计算关联 -----

results2 <- find_associated_metabolites.2(
  phyloseq_obj = ps_aligned,
  metabolome_mat = tab_aligned,
  target_microbe_name = "Arthrobacter",
  top_n = 50
)

names(results2)
results2$top_metabolites

tax_aligned = ps2_aligned %>% tax_table() %>% as.data.frame()
dat_aligned = results2$top_metabolites %>% left_join(tax_aligned, by = c("metabolite"="metab_id"))

p4 <- plot_microbe_circular_network2(dat_aligned,"Microbe_x",thr=0.5)
p4
ggsave(file.path(main_dir, "aligned_cross_omics_network.pdf"), p4, width = 10, height = 10)
ggsave(file.path(main_dir, "aligned_cross_omics_network.png"), p4, width = 10, height = 10, dpi = 300)

# 保存关联代谢物分析结果
associated_metabolites_wb <- createWorkbook()
addWorksheet(associated_metabolites_wb, "associated_metabolites")
addWorksheet(associated_metabolites_wb, "correlation_data")
addWorksheet(associated_metabolites_wb, "top_metabolites_details")
writeData(associated_metabolites_wb, "associated_metabolites", dat_aligned, rowNames = TRUE)
writeData(associated_metabolites_wb, "correlation_data", results2$correlation_data, rowNames = TRUE)
writeData(associated_metabolites_wb, "top_metabolites_details", results2$correlation_results, rowNames = TRUE)
saveWorkbook(associated_metabolites_wb, file.path(main_dir, "associated_metabolites_analysis_results.xlsx"), overwrite = TRUE)

head(dat_aligned)

# smart_metabolite_selection 指定关联代谢物数量 ------

selection_report <- smart_metabolite_selection(
  results = results2,
  sample_size = 20,
  validation_available = FALSE,
  resource_level = "low"
)

# 保存选择报告的图
p_selection1 = selection_report$plots$p1
p_selection2 = selection_report$plots$p2
p_selection_combined = gridExtra::grid.arrange(p_selection1, p_selection2, ncol = 1)

ggsave(file.path(main_dir, "method_comparison.pdf"), plot = p_selection1, width = 12, height = 8, dpi = 300)
ggsave(file.path(main_dir, "method_comparison.png"), plot = p_selection1, width = 12, height = 8, dpi = 300)
ggsave(file.path(main_dir, "score_distribution.pdf"), plot = p_selection2, width = 12, height = 8, dpi = 300)
ggsave(file.path(main_dir, "score_distribution.png"), plot = p_selection2, width = 12, height = 8, dpi = 300)

# 获取推荐的代谢物列表
selected_metabolites <- selection_report$metabolites$standard

dat = selected_metabolites %>% left_join(tax_aligned, by = c("metabolite"="metab_id"))
dat

# 保存智能代谢物选择结果
smart_selection_wb <- createWorkbook()
addWorksheet(smart_selection_wb, "smart_metabolite_selection")
addWorksheet(smart_selection_wb, "selection_report")
addWorksheet(smart_selection_wb, "selection_methods")
writeData(smart_selection_wb, "smart_metabolite_selection", dat, rowNames = TRUE)
writeData(smart_selection_wb, "selection_report", selection_report$report, rowNames = TRUE)
saveWorkbook(smart_selection_wb, file.path(main_dir, "smart_metabolite_selection_results.xlsx"), overwrite = TRUE)

# 创建推荐代谢物的多维评分图 -----

df = dat

# 排序：先正后负；各自内部 integrated_score 降序。让"图中最上面=最重要"
df_ord <- df %>%
  mutate(cor_sign = ifelse(pearson_cor >= 0, "Positive", "Negative"),
         sign_key = ifelse(cor_sign == "Positive", 1L, 0L)) %>%
  arrange(desc(sign_key), desc(integrated_score))
met_levels <- df_ord$Metabolite
df_ord <- df_ord %>%
  mutate(Metabolite = factor(Metabolite, levels = rev(met_levels)))

# 拉成长表（把 pearson_cor 带过去，用于上色）
metric_levels <- c("Integrated","Abundance","Correlation","Network","Differential","ML")
df_long <- df_ord %>%
  dplyr::select(Metabolite, pearson_cor,
         integrated_score, abundance_score, correlation_score,
         network_score, diff_score, ml_score) %>%
  dplyr::rename(
    Integrated   = integrated_score,
    Abundance    = abundance_score,
    Correlation  = correlation_score,
    Network      = network_score,
    Differential = diff_score,
    ML           = ml_score
  ) %>%
  pivot_longer(cols = -c(Metabolite, pearson_cor),
               names_to = "metric", values_to = "score") %>%
  mutate(metric = factor(metric, levels = metric_levels))

# 为了颜色对称，按数据自动取 |r| 的最大值做两端
rmax <- max(abs(df_long$pearson_cor), na.rm = TRUE)

# 左列：Integrated（显示名字；填充=pearson_cor 的红蓝渐变）
p_left_metabolite <- ggplot(dplyr::filter(df_long, metric == "Integrated"),
                            aes(x = score, y = Metabolite, fill = pearson_cor)) +
  geom_col(width = 0.8) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0.01, 0.05))) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-rmax, rmax),
                       name = "Pearson r\n(负→蓝，正→红)") +
  labs(x = "Integrated", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# 右侧多列：其他分数（同样用 pearson_cor 上色；隐藏名字）
p_right_metabolite <- ggplot(dplyr::filter(df_long, metric != "Integrated"),
                             aes(x = score, y = Metabolite, fill = pearson_cor)) +
  geom_col(width = 0.8) +
  facet_grid(cols = vars(metric), scales = "fixed", switch = "x") +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0.01, 0.05))) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-rmax, rmax),
                       name = "Pearson r\n(负→蓝，正→红)") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(face = "bold")
  )

# 拼接，并收集共享图例
p_metabolite_scores <- p_left_metabolite + p_right_metabolite + plot_layout(widths = c(1, 5), guides = "collect") &
  theme(legend.position = "right")


print(p_metabolite_scores)
ggsave(file.path(main_dir, "selected_metabolite_scores.pdf"), p_metabolite_scores, width = 16, height = 10)
ggsave(file.path(main_dir, "selected_metabolite_scores.png"), p_metabolite_scores, width = 16, height = 10, dpi = 300)

# 代谢组学特征挑选 -----

tab_ms_selection = ps.ms %>% subset_samples.wt("Group",c("WT","KO")) %>% vegan_otu() %>%
  as.data.frame()

map = sample_data(ps.ms %>% subset_samples.wt("Group",c("WT","KO")))
map$Group = as.factor(map$Group)
map$ID = NULL

# 运行特征挑选主流程
res_ms <- ms_feature_selection(
  ms_mat         = tab_ms_selection,
  sample_meta    = map,
  group_var      = "Group",
  batch_var      = NULL,
  qc_flag_var    = NULL,
  prevalence_threshold = 0.7,
  impute_method  = "halfmin",
  normalize_method = "median",
  combat         = FALSE,
  qc_drift_correction = TRUE,
  qc_rsd_cut     = 0.30,
  diff_method    = "limma",
  ml_method      = "rf",
  cor_method     = "spearman",
  cor_threshold  = 0.6,
  p_threshold    = 0.05,
  weights = list(
    abundance   = 0.3,
    importance  = 0.3,
    differential= 0.2,
    network     = 0.2
  ),
  nfolds = 5,
  ntree  = 300
)

# 查看结果
names(res_ms)
head(res_ms$final_ranking, 15)

dat_ms = res_ms$final_ranking
tax_ms = ps.ms %>% tax_table() %>% as.data.frame()
dat_ms_final = dat_ms %>% left_join(tax_ms, by = c("Metabolite"="metab_id"))


# 保存代谢组特征选择结果
metabolomics_wb <- createWorkbook()
addWorksheet(metabolomics_wb, "metabolomics_features")
addWorksheet(metabolomics_wb, "selection_statistics")
addWorksheet(metabolomics_wb, "diff_results")
addWorksheet(metabolomics_wb, "importance_results")
addWorksheet(metabolomics_wb, "network_results")
addWorksheet(metabolomics_wb, "abundance_results")
addWorksheet(metabolomics_wb, "processed_data")

writeData(metabolomics_wb, "metabolomics_features", dat_ms_final, rowNames = TRUE)
writeData(metabolomics_wb, "selection_statistics", res_ms$statistics, rowNames = TRUE)
writeData(metabolomics_wb, "diff_results", res_ms$diff_results, rowNames = TRUE)
writeData(metabolomics_wb, "importance_results", res_ms$importance_results, rowNames = TRUE)
writeData(metabolomics_wb, "network_results", res_ms$network_results, rowNames = TRUE)
writeData(metabolomics_wb, "abundance_results", res_ms$abundance_results, rowNames = TRUE)
writeData(metabolomics_wb, "processed_data", res_ms$processed_ms, rowNames = TRUE)
saveWorkbook(metabolomics_wb, file.path(main_dir, "metabolomics_feature_selection_results.xlsx"), overwrite = TRUE)

# 寻找与代谢物相关的微生物 ----

ps2_microbe = ps.ms %>% subset_samples.wt("Group",c("WT","OE"))
tab_microbe = ps2_microbe %>% vegan_otu() %>%
  as.data.frame()

ps_microbe = ps.16s %>% subset_samples.wt("Group",c("WT","OE")) %>% tax_glom_wt(6)

res_microbe <- find_associated_microbes(
  phyloseq_obj = ps_microbe,
  metabolome_mat = tab_microbe,
  target_metabolite_name = "metab_7974",
  top_n = 20,
  cor_method = "spearman",
  prevalence_threshold = 0.10,
  detection_threshold  = 1e-4,
  use_clr = TRUE,
  group_split = "median"
)


# 调整列名以匹配 smart_metabolite_selection 的期望
colnames(res_microbe$correlation_results)[2] = "pearson_cor"
colnames(res_microbe$correlation_results)[4] = "pearson_fdr"

# 保存代谢物-微生物关联结果
metabolite_microbe_wb <- createWorkbook()
addWorksheet(metabolite_microbe_wb, "metabolite_microbe_association")
addWorksheet(metabolite_microbe_wb, "top_microbes_details")
addWorksheet(metabolite_microbe_wb, "correlation_results")
writeData(metabolite_microbe_wb, "metabolite_microbe_association", res_microbe$top_microbes, rowNames = TRUE)
writeData(metabolite_microbe_wb, "top_microbes_details", res_microbe$correlation_matrix, rowNames = TRUE)
writeData(metabolite_microbe_wb, "correlation_results", res_microbe$correlation_results, rowNames = TRUE)
saveWorkbook(metabolite_microbe_wb, file.path(main_dir, "metabolite_microbe_association_results.xlsx"), overwrite = TRUE)

