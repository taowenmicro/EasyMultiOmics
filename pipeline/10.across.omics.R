

library(sva)          # ComBat批次校正
library(RSpectra)     # 稀疏PCA
library(dplyr)        # 数据操作
library(compositions) # CLR
library(ggClusterNet)
library(phyloseq)
library(EasyMultiOmics)
library(caret)
library("DESeq2")
library(limma)
library(randomForest)
library(DESeq2)


#挑选基因组学中的特征#------

ps = ps.16s %>% subset_samples.wt("Group",c("WT","KO"))
# ps = ps %>% subset_taxa.wt("Species",c("unclassified"),TRUE)
ps = ps %>% tax_glom_wt(6)
ps
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

id = dat$OTU[1]
id

#  单一组学特征展示#-------
# ==== 可调参数 ====
top_n <- 50  # 只画前N个（按 Rank 取前N）
df_top <- dat  %>% arrange(Rank) %>% slice_head(n = top_n)
# 必要包
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)


# ==== 排序：先正向(>=0)后负向；组内按 Weighted_score 降序；图中最上面=更重要 ====
df_ord <- df_top %>%
  mutate(sign_key = ifelse(log2FoldChange >= 0, 1L, 0L)) %>%
  arrange(desc(sign_key), desc(Weighted_score))

otu_levels <- df_ord$OTU
df_ord <- df_ord %>% mutate(OTU = factor(OTU, levels = rev(otu_levels)))  # 翻转，使最重要在最上

# ==== 组装长表：左列是 Weighted，右侧面板是其它指标 ====
df_long <- df_ord %>%
  transmute(
    OTU,
    log2FoldChange,
    Weighted       = Weighted_score,
    Mean_abundance = Mean_abundance,
    Prevalence     = Prevalence,
    Significance   = -log10(padj),     # 显著性用 -log10(FDR)
    log2FC         = log2FoldChange,
    Importance     = Importance,
    Degree         = Degree,
    Betweenness    = Betweenness
  ) %>%
  pivot_longer(cols = -c(OTU, log2FoldChange),
               names_to = "metric", values_to = "score")

metric_levels <- c("Weighted","Mean_abundance","Prevalence",
                   "Significance","log2FC","Importance","Degree","Betweenness")
df_long$metric <- factor(df_long$metric, levels = metric_levels)

# ==== 数值标签：按不同列做合适的格式 ====
fmt_sci <- label_scientific(digits = 2)
df_long <- df_long %>%
  mutate(label = dplyr::case_when(
    metric == "Weighted"        ~ number(score, accuracy = 0.001),
    metric == "Mean_abundance"  ~ fmt_sci(score),
    metric == "Prevalence"      ~ percent(score, accuracy = 1),
    metric == "Significance"    ~ number(score, accuracy = 0.1),
    metric == "log2FC"          ~ number(score, accuracy = 0.01),
    metric == "Importance"      ~ number(score, accuracy = 0.1),
    metric == "Degree"          ~ number(score, accuracy = 1),
    metric == "Betweenness"     ~ number(score, accuracy = 0.1),
    TRUE                        ~ as.character(round(score, 3))
  ))

# 颜色范围：按 |log2FC| 的最大值对称
fc_max <- max(abs(df_ord$log2FoldChange), na.rm = TRUE)

# ==== 左列：Weighted（显示OTU名） ====
p_left <- ggplot(dplyr::filter(df_long, metric == "Weighted"),
                 aes(x = score, y = OTU, fill = log2FoldChange)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = label), hjust = -0.1, size = 2.7) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.12))) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-fc_max, fc_max),
                       name = "log2FC\n(蓝负, 红正)") +
  labs(x = "Weighted score", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right"
  )

# ==== 右侧：其余指标（不显示OTU名），每列一个面板，横轴自适应 ====
p_right <- ggplot(dplyr::filter(df_long, metric != "Weighted"),
                  aes(x = score, y = OTU, fill = log2FoldChange)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = label), hjust = -0.1, size = 2.5) +
  facet_grid(cols = vars(metric), scales = "free_x", switch = "x") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.12))) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-fc_max, fc_max),
                       name = "log2FC\n(蓝负, 红正)") +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.placement    = "outside",
    strip.text.x       = element_text(face = "bold"),
    legend.position    = "right"
  )

# ==== 拼接：左1列 + 右多列，共享图例 ====
p <- p_left + p_right + plot_layout(widths = c(1, 6), guides = "collect")
print(p)


# ggsave("../公开发表文章算法验证/干旱SBB验证/微生物特征集合ko.pdf",p,width = 12,height = 8)

# 可选保存
# ggsave("microbe_importance_columns.png", p, width = 18, height = 12, dpi = 300)
#  主图2#-----






# 跨组学特征寻找#-------



ps2 = ps.ms %>% subset_samples.wt("Group",c("WT","KO"))
tab = ps2 %>% vegan_otu() %>%
  as.data.frame()


ps = ps.16s %>% subset_samples.wt("Group",c("WT","KO")) %>% tax_glom_wt(6)

# 运行完整分析
results <- find_associated_metabolites(
  phyloseq_obj = ps,
  metabolome_mat = tab,
  target_microbe_name = "Pseudonocardia" ,
  top_n = 30
)
#
names(results)
results$top_metabolites

tax = ps2 %>% tax_table() %>% as.data.frame()
head(tax)

dat = results$top_metabolites %>% left_join(tax,by = c("metabolite"="metab_id"))

head(dat)


dat$Metabolite
# 查看结果
print(head(results$top_metabolites))

# 跨组学关系寻找可视化#-----

library(dplyr)
library(ggplot2)
library(igraph)
library(tidygraph)
library(ggraph)

# # ===== 输入示例（用你的 dat 替换）=====
# # dat: metabolite, group, r
# # microbe_name: 中心微生物名
# # set.seed(1)
# dat <- tibble::tibble(
#   metabolite = paste0("metab_", 1:60),
#   group = rep(c("G1","G2","G3"), length.out = 60),
#   r = runif(60, -1, 1)
# )


dat = results$top_metabolites %>% left_join(tax,by = c("metabolite"="metab_id"))

head(dat)

id =ifelse(dat$spearman_cor >= 0, 1, -1)

dat1 <- tibble::tibble(
  metabolite =dat$Metabolite,
  group = ifelse(dat$spearman_cor >= 0, "pos", "neg"),
  r = dat$integrated_score*id
)



microbe_name <- "Target_microbe"

thr <- 0.0
dat_use <- dat1 %>% filter(is.finite(r), abs(r) >= thr)

# ===== 组、节点、边（三层：microbe -> group -> metabolite）=====
groups <- dat_use %>%
  distinct(group) %>%
  transmute(name = group, type = "group")

metabs <- dat_use %>%
  arrange(group, desc(abs(r))) %>%
  transmute(name = metabolite, type = "metabolite", group)

microbe <- tibble(name = microbe_name, type = "microbe")

nodes <- bind_rows(microbe, groups, metabs)

edges <- bind_rows(
  groups %>% transmute(from = microbe_name, to = name, r = NA_real_),
  dat_use  %>% transmute(from = group, to = metabolite, r = r)
)

g <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE) %>%
  activate(edges) %>%
  mutate(
    sign  = case_when(is.na(r) ~ "struct", r >= 0 ~ "pos", TRUE ~ "neg"),
    alpha = if_else(sign == "struct", 0.35, 0.65)   # 结构线淡一些
  )
g
# ===== 配色：只区分正/负（结构线灰色）=====
# ---- 仍沿用你上一步构建好的 g（tbl_graph，含 edges: sign/alpha） ----

col_pos <- "#C43E2F"   # 正相关
col_neg <- "#2B6CB0"   # 负相关
col_struct <- "grey80" # 结构线（microbe->group）

label_offset <- 0.05   # 标签沿半径外移距离（可调）

p <- ggraph(g, layout = "dendrogram", circular = TRUE) +
  # 细、仅按正负上色；结构线淡灰
  geom_edge_diagonal(
    aes(edge_colour = sign, edge_alpha = alpha),
    edge_width = 0.28, lineend = "round"
  ) +
  scale_edge_colour_manual(
    values = c(pos = col_pos, neg = col_neg, struct = col_struct),
    breaks = c("pos","neg"),
    labels = c("Positive", "Negative"),
    name   = "Correlation"
  ) +
  guides(edge_alpha = "none") +

  # 节点（不改）
  geom_node_point(
    aes(shape = type, fill = type),
    size = 2, colour = "white", stroke = 0.5
  ) +
  scale_shape_manual(values = c(microbe = 21, group = 21, metabolite = 21), guide = "none") +
  scale_fill_manual(values = c(microbe = "#8C564B", group = "#EADBC8", metabolite = "#C43E2F"), guide = "none") +

  # 代谢物标签：先算偏移后坐标与角度，再映射到 x/y/angle/hjust
  geom_node_text(
    data = function(n){
      n %>%
        dplyr::filter(type == "metabolite") %>%
        dplyr::mutate(
          ang    = atan2(y, x) * 180 / pi,
          flip   = ang > 90 & ang < 270,
          lab_ang = ifelse(flip, ang - 180, ang),
          lab_h   = ifelse(flip, 1, 0),
          r     = sqrt(x^2 + y^2),
          ux    = ifelse(r > 0, x/r, 0),
          uy    = ifelse(r > 0, y/r, 0),
          x_lab = x + ux * label_offset,
          y_lab = y + uy * label_offset
        )
    },
    aes(x = x_lab, y = y_lab, label = name, angle = lab_ang, hjust = lab_h),
    size = 2.2, family = "sans", check_overlap = TRUE
  ) +

  # 组名与中心标签（不改）
  geom_node_label(
    data = function(n) dplyr::filter(n, type == "group"),
    aes(label = name), label.size = 0, size = 3.0, fill = "white",alpha = 0
  ) +
  geom_node_label(
    data = function(n) dplyr::filter(n, type == "microbe"),
    aes(label = name), label.size = 0, size = 3.6, fill = "white",alpha = 0
  ) +

  coord_equal(clip = "off") +
  theme_void(base_size = 11) +
  theme(
    legend.position = "right",
    plot.margin = grid::unit(c(60, 60, 60, 60), "pt")
  )

p
# ggsave("./cs.pdf",p,width = 5,height = 6)


#  multi_omics_alignment 多组学矫正#-----



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
# 基本使用
alignment_results <- multi_omics_alignment(
  microbiome_data = otu,  # 样本×微生物
  metabolome_data = tab,  # 样本×代谢物
  n_factors = 10,                       # 潜在因子数
  method = "fast"                       # 自动选择快速或完整方法
)

# 查看对齐质量
print(alignment_results$quality_details)
names(alignment_results)
alignment_results$aligned_microbiome
alignment_results$aligned_metabolome

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
library(igraph)
# find_associated_metabolites.2 计算关联#-----
# 运行完整分析
results <- find_associated_metabolites.2(
  phyloseq_obj = ps,
  metabolome_mat = tab,
  target_microbe_name = "Pseudonocardia",
  top_n = 50
)

names(results)
results$top_metabolites

tax = ps2 %>% tax_table() %>% as.data.frame()
head(tax)

dat = results$top_metabolites %>% left_join(tax,by = c("metabolite"="metab_id"))

head(dat)
dat %>% tail()
dat$Metabolite

results$correlation_results
# smart_metabolite_selection 指定关联代谢物数量#------
selection_report <- smart_metabolite_selection(
  results = results,
  sample_size = 20,           # 你的样本量
  validation_available = F, # 是否有验证队列
  resource_level = "low"    # 资源水平: "low", "medium", "high"
)

# 获取推荐的代谢物列表
selected_metabolites <- selection_report$metabolites$standard

tax = ps2 %>% tax_table() %>% as.data.frame()
head(tax)

dat = selected_metabolites %>% left_join(tax,by = c("metabolite"="metab_id"))
dat$Metabolite
head(dat)



p = gridExtra::grid.arrange(selection_report$plots$p1,selection_report$plots$p2, ncol = 1)# 需要的包

# ggsave("./cs.pdf",p,width = 9,height = 7)




library(ggplot2)
library(dplyr)
library(tidyr)

# 假设你的数据框名为 df（即你贴出来的那张表）
# df <- read.csv("your_table.csv")
df =selected_metabolites


tax = ps2 %>% tax_table() %>% as.data.frame()
head(tax)

dat = selected_metabolites %>% left_join(tax,by = c("metabolite"="metab_id"))
dat$Metabolite
head(dat)
df = dat



library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# 1) 排序：先正后负；各自内部 integrated_score 降序。让“图中最上面=最重要”
df_ord <- df %>%
  mutate(cor_sign = ifelse(pearson_cor >= 0, "Positive", "Negative"),
         sign_key = ifelse(cor_sign == "Positive", 1L, 0L)) %>%
  arrange(desc(sign_key), desc(integrated_score))
met_levels <- df_ord$Metabolite
df_ord <- df_ord %>%
  mutate(Metabolite = factor(Metabolite, levels = rev(met_levels)))  # 最高在最上面

# 2) 拉成长表（把 pearson_cor 带过去，用于上色）
metric_levels <- c("Integrated","Abundance","Correlation","Network","Differential","ML")
df_long <- df_ord %>%
  select(Metabolite, pearson_cor,
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

# 3) 左列：Integrated（显示名字；填充=pearson_cor 的红蓝渐变）
p_left <- ggplot(dplyr::filter(df_long, metric == "Integrated"),
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

# 4) 右侧多列：其他分数（同样用 pearson_cor 上色；隐藏名字）
p_right <- ggplot(dplyr::filter(df_long, metric != "Integrated"),
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

# 5) 拼接，并收集共享图例
p <- p_left + p_right + plot_layout(widths = c(1, 5), guides = "collect") &
  theme(legend.position = "right")

print(p)

# ggsave("./cs2.pdf",p,width = 12,height = 7)

# ggsave("top50_sign_gradient_columns.png", p, width = 16, height = 10, dpi = 300)


#  代谢组学特征挑选
## =========================
tab = ps.ms %>%subset_samples.wt("Group",c("WT","KO")) %>% vegan_otu( ) %>%
  as.data.frame()

map = sample_data(ps.ms %>%subset_samples.wt("Group",c("WT","KO")) )
map$Group = as.factor(map$Group)
map$ID = NULL
## 2) 运行特征挑选主流程
res <- ms_feature_selection(
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
# head(res$diff_results)          # limma 结果
# head(res$importance_results)    # ML 重要性
# head(res$network_results)       # 网络指标
# head(res$abundance_results)     # 丰度/流行率
# dim(res$processed_ms)           # 处理后矩阵（可用于下游分析）


dat = res$final_ranking
head(dat)

tax = ps.ms %>% tax_table() %>% as.data.frame()
head(tax)

dat =  res$final_ranking %>% left_join(tax,by = c("Metabolite"="metab_id"))

head(dat)
dat$Metabolite.y %>% head(50)

dat$Metabolite.y[1]





ps2 = ps.ms %>% subset_samples.wt("Group",c("WT","OE"))
tab = ps2 %>% vegan_otu() %>%
  as.data.frame()


ps = ps.16s %>% subset_samples.wt("Group",c("WT","OE")) %>% tax_glom_wt(6)


# phyloseq_obj = ps
# metabolome_mat = tab
# target_metabolite_name = "metab_4814"  # 我们植入信号的代谢物
# top_n = 20
# cor_method = "spearman"
# prevalence_threshold = 0.10
# detection_threshold  = 1e-4
# use_clr = TRUE
# group_split = "median"


res_sim <- find_associated_microbes(
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
res_sim$top_microbes

colnames(res_sim$correlation_results)[2] = "pearson_cor"
colnames(res_sim$correlation_results)[4] = "pearson_fdr"



# smart_metabolite_selection 指定关联代谢物数量#------
selection_report <- smart_metabolite_selection(
  results = res_sim,
  sample_size = 20,           # 你的样本量
  validation_available = F, # 是否有验证队列
  resource_level = "low"    # 资源水平: "low", "medium", "high"
)

# 获取推荐的代谢物列表
selected_metabolites <- selection_report$metabolites$standard



