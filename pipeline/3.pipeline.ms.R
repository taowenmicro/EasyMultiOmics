# 代谢组学数据挖掘-----
rm(list=ls())
library(EasyMultiOmics)
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(ggrepel)
library(EasyMultiOmics.db)
library(openxlsx)

## ===================== 0. 基础设置 & 通用函数 =====================


## ===================== 1. 读入数据 & 主题设置 =====================
# 方式一：交互选择
# file_path <- tcltk::tk_choose.files(
#   caption = "请选择 ps_ITS.rds 文件",
#   multi   = FALSE,
#   filters = matrix(c("RDS files", ".rds",
#                      "All files", "*"), ncol = 2, byrow = TRUE)
# )
# 方式二：直接指定路径
# file_path <- "./data/ps_16s.rds"
# ps.16s    <- readRDS(file_path)
# ps.16s
# ps.ms <- EasyMultiOmics::ps.ms
ps.ms <-readRDS( "E:/Shared_Folder/Help_project24.01/冯明峰病毒番茄微生物/data/ps_GC.rds") %>% subset_taxa.wt("KEGG COMPOUND ID",NA,TRUE)


tax = ps.ms %>% vegan_tax() %>% as.data.frame()
head(tax)

# # 基础路径
# metabolite_path <- "../result/metabolite/"
# dir.create(metabolite_path, recursive = TRUE, showWarnings = FALSE)
out_dir <- create_omics_result_dir_auto(ps.ms,
                                        base_dir = "../result",
                                        include_time = FALSE)



gnum <- phyloseq::sample_data(ps.ms)$Group %>% unique() %>% length()
axis_order <- phyloseq::sample_data(ps.ms)$Group %>% unique()
col.g <- get_group_cols(axis_order)
scales::show_col(col.g)



package.amp()
res <- theme_my(ps.ms)
mytheme1 <- res[[1]]
mytheme2 <- res[[2]]
colset1 <- res[[3]]
colset2 <- res[[4]]
colset3 <- res[[5]]
colset4 <- res[[6]]

## ===================== 2. 注释 & 标准化 =====================

metabolite_annotation_path <- file.path(metabolite_path, "01_annotation")
dir.create(metabolite_annotation_path, recursive = TRUE, showWarnings = FALSE)

annotation_xlsx_path <- file.path(metabolite_annotation_path, "annotation_results.xlsx")

if (file.exists(annotation_xlsx_path)) {
  metabolite_annotation_wb <- openxlsx::loadWorkbook(annotation_xlsx_path)
} else {
  metabolite_annotation_wb <- openxlsx::createWorkbook()
}

### ---------- 1. ann.HMDB: 代谢物注释 ----------

tax <- ps.ms %>% vegan_tax() %>%
  as.data.frame()

head(tax)
id <- tax$ID

tax1 <- ann.HMDB(id = id)


tax1 <- tax1 %>% distinct(id, .keep_all = TRUE) %>%
  column_to_rownames("id")

tax$ID <- row.names(tax)
tax3 <- tax %>% left_join(tax1, by = c("ID" = "id.org"))
row.names(tax3) <- tax3$metab_id
tax_table(ps.ms) <- phyloseq::tax_table(as.matrix(tax3))

write_sheet2(metabolite_annotation_wb, "HMDB_annotation", tax3)
openxlsx::saveWorkbook(metabolite_annotation_wb, annotation_xlsx_path, overwrite = TRUE)

### ---------- 2. ann.kegg2: 代谢物注释KEGG数据库 ----------

tax2 <- ann.kegg(id)

write_sheet2(metabolite_annotation_wb, "KEGG_annotation", tax2)
openxlsx::saveWorkbook(metabolite_annotation_wb, annotation_xlsx_path, overwrite = TRUE)

### ---------- 3. 数据预处理 ----------

ps.ms2 <- zone.fill.ms(ps = ps.ms, method = "repeat")
# ps.ms2 <- scale_IS.ms(ps = ps.ms, IS = "metab_8497")
# ps.ms2 <- scale_QC.ms(ps = ps.ms, QC = c("OE36_2", "OE36_3", "OE2_1"))
ps.ms2 <- normalize.ms(ps = ps.ms, method = "rela")

## ===================== 3. 排序分析 =====================

metabolite_ordinate_path <- file.path(metabolite_path, "02_ordinate")
dir.create(metabolite_ordinate_path, recursive = TRUE, showWarnings = FALSE)

ordinate_xlsx_path <- file.path(metabolite_ordinate_path, "ordinate_results.xlsx")

if (file.exists(ordinate_xlsx_path)) {
  metabolite_ordinate_wb <- openxlsx::loadWorkbook(ordinate_xlsx_path)
} else {
  metabolite_ordinate_wb <- openxlsx::createWorkbook()
}

### ---------- 7. ordinate.ms: 代谢物排序分析 ----------

result <- ordinate.ms(ps = ps.ms, group = "Group", dist = "bray",
                      method = "PCoA", Micromet = "anosim", pvalue.cutoff = 0.05,
                      pair = FALSE)

p3_1 <- result[[1]] +
  scale_fill_manual(values = colset1) +
  scale_color_manual(values = colset1, guide = FALSE) +
  mytheme1

plotdata <- result[[2]]

p3_2 <- result[[3]] +
  scale_fill_manual(values = colset1) +
  scale_color_manual(values = colset1, guide = FALSE) +
  mytheme1

cent <- aggregate(cbind(x, y) ~ Group, data = plotdata, FUN = mean)
segs <- merge(plotdata, setNames(cent, c('Group', 'oNMDS1', 'oNMDS2')),
              by = 'Group', sort = FALSE)

p3_3 <- p3_1 + geom_segment(data = segs,
                            mapping = aes(xend = oNMDS1,
                                          yend = oNMDS2, color = Group), show.legend = FALSE,
                            alpha = 0.5
) +
  geom_point(mapping = aes(x = x, y = y, fill = Group),
             data = cent, size = 5, pch = 23,
             color = "black") +
  scale_fill_manual(values = colset1) +
  scale_color_manual(values = colset1, guide = FALSE) +
  mytheme1

save_plot2(p3_1, metabolite_ordinate_path, "PCoA_basic", width = 8, height = 6)
save_plot2(p3_2, metabolite_ordinate_path, "PCoA_labeled", width = 8, height = 6)
save_plot2(p3_3, metabolite_ordinate_path, "PCoA_refined", width = 8, height = 6)

write_sheet2(metabolite_ordinate_wb, "PCoA_coordinates", plotdata)
write_sheet2(metabolite_ordinate_wb, "centroids", cent)
openxlsx::saveWorkbook(metabolite_ordinate_wb, ordinate_xlsx_path, overwrite = TRUE)

### ---------- 8. MicroTest.ms: 代谢组总体差异检测 ----------

dat1 <- MicroTest.ms(ps = ps.ms, Micromet = "adonis", dist = "bray")

write_sheet2(metabolite_ordinate_wb, "overall_test", dat1)
openxlsx::saveWorkbook(metabolite_ordinate_wb, ordinate_xlsx_path, overwrite = TRUE)

### ---------- 9. pairMicroTest.ms: 两两分组代谢总体水平差异检测 ----------

dat2 <- pairMicroTest.ms(ps = ps.ms, Micromet = "MRPP", dist = "bray")

write_sheet2(metabolite_ordinate_wb, "pairwise_test", dat2)
openxlsx::saveWorkbook(metabolite_ordinate_wb, ordinate_xlsx_path, overwrite = TRUE)

### ---------- 10. mantal.ms: 代谢群落差异检测普鲁士分析 ----------

result <- mantal.ms(ps = ps.ms,
                    method = "spearman",
                    group = "Group",
                    ncol = gnum,
                    nrow = 1)

data <- result[[1]]
p3_7 <- result[[2]] + mytheme1

save_plot2(p3_7, metabolite_ordinate_path, "mantel_test", width = gnum*5, height = 6)

write_sheet2(metabolite_ordinate_wb, "mantel_test", data)
openxlsx::saveWorkbook(metabolite_ordinate_wb, ordinate_xlsx_path, overwrite = TRUE)

### ---------- 11. plsda.ms: plsda排序分析 ----------

library(ggalt)
library(vegan)
library(mixOmics)
library(ggforce)
library(caret)

res <- plsda_ms(ps = ps.ms, Group = "Group")
p11 <- res[[1]]
dat <- res[[2]]

save_plot2(p11, metabolite_ordinate_path, "PLSDA_plot", width = 8, height = 6)

write_sheet2(metabolite_ordinate_wb, "PLSDA_data", dat)
openxlsx::saveWorkbook(metabolite_ordinate_wb, ordinate_xlsx_path, overwrite = TRUE)

### ---------- 12. oplsda.ms: oplsda排序分析 ----------

library(ropls)

res <- oplsda.ms(ps = ps.ms, ncol = 3, nrow = 1)

p12 <- res[1]
dat <- res[2]

if(is.ggplot(p12[[1]])) {
  save_plot2(p12[[1]], metabolite_ordinate_path, "OPLSDA_plot", width = 12, height = 4)
}
if(is.data.frame(dat[[1]])) {
  write_sheet2(metabolite_ordinate_wb, "OPLSDA_data", dat[[1]])
}
openxlsx::saveWorkbook(metabolite_ordinate_wb, ordinate_xlsx_path, overwrite = TRUE)

## ===================== 4. 代谢物分类 =====================

metabolite_classification_path <- file.path(metabolite_path, "03_classification")
dir.create(metabolite_classification_path, recursive = TRUE, showWarnings = FALSE)

classification_xlsx_path <- file.path(metabolite_classification_path, "classification_results.xlsx")

if (file.exists(classification_xlsx_path)) {
  metabolite_classification_wb <- openxlsx::loadWorkbook(classification_xlsx_path)
} else {
  metabolite_classification_wb <- openxlsx::createWorkbook()
}

### ---------- 13. Ven.Upset.ms: 用于展示共有、特有的代谢物 ----------

res <- Ven.Upset.metm(ps = ps.ms,
                      group = "Group",
                      N = 0.5,
                      size = 3)
p10.1 <- res[[1]]
dat <- res[[3]]
head(dat)
str(p10.1)

save_plot2(p10.1, metabolite_classification_path, "Venn_Upset_plot", width = 10, height = 8)

write_sheet2(metabolite_classification_wb, "Venn_Upset_data", dat)
openxlsx::saveWorkbook(metabolite_classification_wb, classification_xlsx_path, overwrite = TRUE)

### ---------- 14. ggflower.ms: 花瓣图展示共有特有代谢物 ----------

res <- ggflower.ms(ps = ps.ms,
                   group = "Group",
                   start = 1,
                   m1 = 1,
                   a = 0.2,
                   b = 1,
                   lab.leaf = 1,
                   col.cir = "yellow",
                   N = 0.5)

p14 <- res[[1]]
dat <- res[[2]]

save_plot2(p14, metabolite_classification_path, "flower_plot", width = 8, height = 8)

write_sheet2(metabolite_classification_wb, "flower_data", dat)
openxlsx::saveWorkbook(metabolite_classification_wb, classification_xlsx_path, overwrite = TRUE)

### ---------- 15. Ms_tern.ms: 三元图展示组成 ----------

ps1 <- ps.ms %>% subset_samples.wt("Group","GM",TRUE) %>%
  filter_OTU_ps(500)
res <- Ms_tern.ms(ps1)
p15 <- res[[1]]
p15_final <- p15[[1]] + theme_bw()

dat <- res[[2]]

save_plot2(p15_final, metabolite_classification_path, "ternary_plot", width = 8, height = 8)

write_sheet2(metabolite_classification_wb, "ternary_data", dat)
openxlsx::saveWorkbook(metabolite_classification_wb, classification_xlsx_path, overwrite = TRUE)

### ---------- 16. barMainplot.ms: 代谢物分类堆叠柱状图 ----------

j <- "Class"
result <- barMainplot.ms(ps = ps.ms,
                         j = "Class",
                         label = FALSE,
                         sd = FALSE,
                         Top = 12)
p4_1 <- result[[1]] + scale_fill_brewer(palette = "Paired")
p4_2 <- result[[3]] +
  scale_fill_manual(values = colset3) +
  scale_x_discrete(limits = axis_order)

databar <- result[[2]] %>%
  dplyr::group_by(Group, aa) %>%
  dplyr::summarise(sum(Abundance)) %>% as.data.frame()
colnames(databar) <- c("Group", j, "Abundance(%)")

save_plot2(p4_1, metabolite_classification_path, "stacked_barplot1", width = 10, height = 6)
save_plot2(p4_2, metabolite_classification_path, "stacked_barplot2", width = 10, height = 6)

write_sheet2(metabolite_classification_wb, "stacked_bar_data", databar)
openxlsx::saveWorkbook(metabolite_classification_wb, classification_xlsx_path, overwrite = TRUE)

### ---------- 17. cluMicro.bar.ms: 聚类堆积柱状图展示组成 ----------

res <- cluMicro.bar.ms(dist = "bray",
                       ps = ps.ms,
                       j = "Class",
                       Top = 10,
                       tran = TRUE,
                       hcluter_method = "complete",
                       Group = "Group",
                       cuttree = length(unique(phyloseq::sample_data(ps.ms)$Group)))

p17.1 <- res[[1]]
p17.2 <- res[[2]]
p17.3 <- res[[3]]
p17.4 <- res[[4]]
clubardata <- res[[5]]

save_plot2(p17.1, metabolite_classification_path, "cluster_barplot1", width = 10, height = 6)
save_plot2(p17.2, metabolite_classification_path, "cluster_barplot2", width = 10, height = 6)
save_plot2(p17.3, metabolite_classification_path, "cluster_barplot3", width = 10, height = 6)
save_plot2(p17.4, metabolite_classification_path, "cluster_barplot4", width = 10, height = 6)

write_sheet2(metabolite_classification_wb, "cluster_bar_data", clubardata)
openxlsx::saveWorkbook(metabolite_classification_wb, classification_xlsx_path, overwrite = TRUE)

## ===================== 5. 差异分析 =====================

metabolite_difference_path <- file.path(metabolite_path, "04_difference")
dir.create(metabolite_difference_path, recursive = TRUE, showWarnings = FALSE)

difference_xlsx_path <- file.path(metabolite_difference_path, "difference_results.xlsx")

if (file.exists(difference_xlsx_path)) {
  metabolite_difference_wb <- openxlsx::loadWorkbook(difference_xlsx_path)
} else {
  metabolite_difference_wb <- openxlsx::createWorkbook()
}

### ---------- 18. cluster_plot.ms: 代谢物层次聚类 ----------

res <- cluster_plot.ms(ps = ps.ms,
                       hcluter_method = "complete",
                       dist = "bray", cuttree = gnum,
                       row_cluster = TRUE,
                       col_cluster = TRUE)

p0 <- res[[1]]
p1 <- res[[2]]
p2 <- res[[3]]
dat <- res[4]

save_plot2(p0, metabolite_difference_path, "cluster_plot1", width = 10, height = 8)
save_plot2(p1, metabolite_difference_path, "cluster_plot2", width = 10, height = 8)
save_plot2(p2, metabolite_difference_path, "cluster_plot3", width = 10, height = 8)

write_sheet2(metabolite_difference_wb, "cluster_data", dat[[1]])
openxlsx::saveWorkbook(metabolite_difference_wb, difference_xlsx_path, overwrite = TRUE)

### ---------- 19. heatmap.ms: 热图展示代谢物差异 ----------

ps.ms_rela <- ps.ms %>% scale_micro(method = "rela") %>%
  tax_glom_wt(ranks = "Class")
result <- heatmap.ms(ps_rela = ps.ms_rela,
                     label = TRUE,
                     col_cluster = TRUE,
                     row_cluster = TRUE)
p19 <- result[[1]]
p19.1 <- result[[2]]
dat <- result[[3]]

save_plot2(p19, metabolite_difference_path, "heatmap1", width = 10, height = 8)
save_plot2(p19.1, metabolite_difference_path, "heatmap2", width = 10, height = 8)

write_sheet2(metabolite_difference_wb, "heatmap_data", dat)
openxlsx::saveWorkbook(metabolite_difference_wb, difference_xlsx_path, overwrite = TRUE)

### ---------- 20. statSuper: 差异代谢物 ----------

result1 <- statSuper(ps = ps.ms, group = "Group", artGroup = NULL, method = "wilcox")
result2 <- statSuper(ps = ps.ms, group = "Group", artGroup = NULL, method = "ttext")

write_sheet2(metabolite_difference_wb, "wilcox_test", result1)
write_sheet2(metabolite_difference_wb, "t_test", result2)
openxlsx::saveWorkbook(metabolite_difference_wb, difference_xlsx_path, overwrite = TRUE)

### ---------- 21. MuiKwWlx2: 分类化合物分组差异 ----------

library(EasyStat)

map <- sample_data(ps.ms)
map <- map[, 1:2]
sample_data(ps.ms) <- map

dat <- ps.ms %>% scale_micro(method = "rela") %>%
  tax_glom_wt(ranks = "Class") %>%
  vegan_otu() %>%
  as.data.frame()

dat$id <- row.names(dat)

dat2 <- dat %>%
  dplyr::left_join(as.tibble(sample_data(ps.ms)), by = c("id" = "ID")) %>%
  dplyr::rename(group = Group) %>%
  dplyr::select(id, group, everything())

dat2$group <- as.factor(dat2$group)

result <- MuiKwWlx2(data = dat2, num = c(3:dim(dat2)[2]))

result1 <- FacetMuiPlotresultBox(data = dat2, num = c(3:dim(dat2)[2]), result = result, sig_show = "abc",
                                 ncol = 4)
p1_1 <- result1[[1]] +
  mytheme2 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

res <- FacetMuiPlotresultBar(data = dat2, num = c(3:dim(dat2)[2]), result = result, sig_show = "abc",
                             ncol = 4)
p1_2 <- res[[1]] +
  guides(color = FALSE) +
  mytheme2 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

res <- FacetMuiPlotReBoxBar(data = dat2, num = c(3:dim(dat2)[2]), result = result, sig_show = "abc", ncol = 4)
p1_3 <- res[[1]] +
  mytheme2 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

save_plot2(p1_1, metabolite_difference_path, "class_boxplot", width = 12, height = 8)
save_plot2(p1_2, metabolite_difference_path, "class_barplot", width = 12, height = 8)
save_plot2(p1_3, metabolite_difference_path, "class_boxbar", width = 12, height = 8)

write_sheet2(metabolite_difference_wb, "class_difference", result)
openxlsx::saveWorkbook(metabolite_difference_wb, difference_xlsx_path, overwrite = TRUE)

### ---------- 22. FacetMuiPlotresultBox: 单变量统计分析-箱线图可视化 ----------

dat <- ps.ms %>%
  ggClusterNet::vegan_otu() %>%
  as.data.frame()

map <- sample_data(ps.ms) %>% as.tibble() %>%
  dplyr::select(ID, Group)
data <- cbind(map[, c(1, 2)], dat)
colnames(data)[2] <- "group"

num <- c(3:ncol(data))

n.fac <- length(num) / 25
n.fac2 <- ceiling(n.fac)
A <- list()

for (j in 1:n.fac2) {
  if (j == 1) {
    A[[j]] <- num[1:25]
  } else if(j != n.fac2){
    x <- (25*(j - 1) + 1)
    y <- 25*j
    A[[j]] <- num[x:y]
  }else if (j == n.fac2){
    x <- (25*(j - 1) + 1)
    y <- 25*j
    A[[j]] <- num[x:length(num)]
  }
}

for (i in 1:n.fac2) {
  result <- EasyStat::MuiaovMcomper2(data = data, num = A[[i]])

  result1 <- EasyStat::FacetMuiPlotresultBox(data = data, num = A[[i]],
                                             result = result,
                                             sig_show = "abc", ncol = 5)
  p1_1 <- result1[[1]] +
    ggplot2::scale_x_discrete(limits = axis_order) +
    mytheme2 +
    ggplot2::guides(fill = guide_legend(title = NULL)) +
    ggplot2::scale_fill_manual(values = colset1)

  save_plot2(p1_1, metabolite_difference_path, paste0("univariate_boxplot_batch", i), width = 15, height = 10)

  write_sheet2(metabolite_difference_wb, paste0("univariate_stats_batch", i), result)
}

openxlsx::saveWorkbook(metabolite_difference_wb, difference_xlsx_path, overwrite = TRUE)

### ---------- 23. FacetMuiPlotresultBar: 单变量统计分析-柱状图可视化 ----------

dat <- ps.ms %>%
  ggClusterNet::vegan_otu() %>%
  as.data.frame()

map <- sample_data(ps.ms) %>% as.tibble() %>%
  dplyr::select(ID, Group)
data <- cbind(map[, c(1, 2)], dat)
colnames(data)[2] <- "group"

num <- c(3:ncol(data))

n.fac <- length(num) / 25
n.fac2 <- ceiling(n.fac)
A <- list()

for (j in 1:n.fac2) {
  if (j == 1) {
    A[[j]] <- num[1:25]
  } else if(j != n.fac2){
    x <- (25*(j - 1) + 1)
    y <- 25*j
    A[[j]] <- num[x:y]
  }else if (j == n.fac2){
    x <- (25*(j - 1) + 1)
    y <- 25*j
    A[[j]] <- num[x:length(num)]
  }
}

for (i in 1:n.fac2) {
  result <- EasyStat::MuiaovMcomper2(data = data, num = A[[i]])

  res <- EasyStat::FacetMuiPlotresultBar(data = data, num = A[[i]],
                                         result = result, sig_show = "abc", ncol = 5)

  p1_2 <- res[[1]] +
    scale_x_discrete(limits = axis_order) +
    guides(color = FALSE) +
    mytheme2 +
    guides(fill = guide_legend(title = NULL)) +
    scale_fill_manual(values = colset1)

  save_plot2(p1_2, metabolite_difference_path, paste0("univariate_barplot_batch", i), width = 15, height = 10)
}

### ---------- 24. FacetMuiPlotReBoxBar: 单变量统计分析-柱状图结合散点图可视化 ----------

dat <- ps.ms %>%
  ggClusterNet::vegan_otu() %>%
  as.data.frame()

map <- sample_data(ps.ms) %>% as.tibble() %>%
  dplyr::select(ID, Group)
data <- cbind(map[, c(1, 2)], dat)
colnames(data)[2] <- "group"

num <- c(3:ncol(data))

n.fac <- length(num) / 25
n.fac2 <- ceiling(n.fac)
A <- list()

for (j in 1:n.fac2) {
  if (j == 1) {
    A[[j]] <- num[1:25]
  } else if(j != n.fac2){
    x <- (25*(j - 1) + 1)
    y <- 25*j
    A[[j]] <- num[x:y]
  }else if (j == n.fac2){
    x <- (25*(j - 1) + 1)
    y <- 25*j
    A[[j]] <- num[x:length(num)]
  }
}

for (i in 1:n.fac2) {
  result <- EasyStat::MuiaovMcomper2(data = data, num = A[[i]])

  res <- EasyStat::FacetMuiPlotReBoxBar(data = data, num = A[[i]],
                                        result = result, sig_show = "abc", ncol = 5)
  p1_3 <- res[[1]] +
    scale_x_discrete(limits = axis_order) +
    mytheme2 +
    guides(fill = guide_legend(title = NULL)) +
    scale_fill_manual(values = colset1)

  save_plot2(p1_3, metabolite_difference_path, paste0("univariate_boxbar_batch", i), width = 15, height = 10)
}

### ---------- 25. MuiHeatmapBubplot: 单变量统计分析-气泡图可视化 ----------

dat <- ps.ms %>%
  ggClusterNet::vegan_otu() %>%
  as.data.frame()

map <- sample_data(ps.ms) %>% as.tibble() %>%
  dplyr::select(ID, Group)
data <- cbind(map[, c(1, 2)], dat)
colnames(data)[2] <- "group"

num <- c(3:ncol(data))

n.fac <- length(num) / 25
n.fac2 <- ceiling(n.fac)
A <- list()

for (j in 1:n.fac2) {
  if (j == 1) {
    A[[j]] <- num[1:25]
  } else if(j != n.fac2){
    x <- (25*(j - 1) + 1)
    y <- 25*j
    A[[j]] <- num[x:y]
  }else if (j == n.fac2){
    x <- (25*(j - 1) + 1)
    y <- 25*j
    A[[j]] <- num[x:length(num)]
  }
}

for (i in 1:n.fac2) {
  result <- EasyStat::MuiaovMcomper2(data = data, num = A[[i]])

  res <- EasyStat::MuiHeatmapBubplot(
    data = data,
    i = A[[i]],
    col_cluster = FALSE,
    row_cluster = FALSE,
    label = TRUE,
    result = result,
    sample = TRUE,
    scale = TRUE
  )
  p1 <- res[[1]]
  p2 <- res[[2]]

  save_plot2(p1, metabolite_difference_path, paste0("heatmap_sample_batch", i), width = 12, height = 8)
  save_plot2(p2, metabolite_difference_path, paste0("bubble_sample_batch", i), width = 12, height = 8)

  res <- EasyStat::MuiHeatmapBubplot(
    data = data,
    i = A[[i]],
    result = result,
    col_cluster = FALSE,
    row_cluster = FALSE,
    label = TRUE,
    sample = FALSE,
    scale = TRUE
  )

  p1 <- res[[1]]
  p2 <- res[[2]]

  save_plot2(p1, metabolite_difference_path, paste0("heatmap_group_batch", i), width = 12, height = 8)
  save_plot2(p2, metabolite_difference_path, paste0("bubble_group_batch", i), width = 12, height = 8)
}

### ---------- 26. value_stackBar: 单变量统计分析-堆叠柱状图可视化 ----------

dat <- ps.ms %>%
  ggClusterNet::vegan_otu() %>%
  as.data.frame()

map <- sample_data(ps.ms) %>% as.tibble() %>%
  dplyr::select(ID, Group)
data <- cbind(map[, c(1, 2)], dat)
colnames(data)[2] <- "group"

num <- c(3:ncol(data))

n.fac <- length(num) / 25
n.fac2 <- ceiling(n.fac)
A <- list()

for (j in 1:n.fac2) {
  if (j == 1) {
    A[[j]] <- num[1:25]
  } else if(j != n.fac2){
    x <- (25*(j - 1) + 1)
    y <- 25*j
    A[[j]] <- num[x:y]
  }else if (j == n.fac2){
    x <- (25*(j - 1) + 1)
    y <- 25*j
    A[[j]] <- num[x:length(num)]
  }
}

for (i in 1:n.fac2) {
  result <- EasyStat::MuiaovMcomper2(data = data, num = A[[i]])

  res <- EasyStat::value_stackBar(
    data = data,
    i = A[[i]],
    result = result,
    add_abc = TRUE)

  p1 <- res[[1]]

  save_plot2(p1, metabolite_difference_path, paste0("stack_bar_batch", i), width = 12, height = 8)
}
## ===================== 6. 生物标志物识别 =====================

metabolite_biomarker_path <- file.path(metabolite_path, "05_biomarker")
dir.create(metabolite_biomarker_path, recursive = TRUE, showWarnings = FALSE)

biomarker_xlsx_path <- file.path(metabolite_biomarker_path, "biomarker_results.xlsx")

if (file.exists(biomarker_xlsx_path)) {
  metabolite_biomarker_wb <- openxlsx::loadWorkbook(biomarker_xlsx_path)
} else {
  metabolite_biomarker_wb <- openxlsx::createWorkbook()
}

id <- sample_data(ps.ms)$Group %>% unique()
aaa <- combn(id, 2)
i <- 1
group <- c(aaa[1, i], aaa[2, i])

pst <- ps.ms %>% subset_samples.wt("Group", group) %>%
  filter_taxa(function(x) sum(x) > 10, TRUE)

### ---------- 27. rfcv.ms: 交叉验证结果 ----------

library(randomForest)
library(caret)
library(ROCR)
library(e1071)

result <- rfcv.ms(ps = ps.ms, group = "Group", optimal = 20, nrfcvnum = 6)
prfcv <- result[[1]]
rfcvtable <- result[[3]]

save_plot2(prfcv, metabolite_biomarker_path, "rfcv_plot", width = 8, height = 6)

write_sheet2(metabolite_biomarker_wb, "rfcv_results", rfcvtable)
openxlsx::saveWorkbook(metabolite_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 28. randomforest.ms: 筛选特征代谢物 ----------

res <- randomforest.ms(ps = ps.ms, group = "Group", optimal = 50)
p25 <- res[[1]]
dat <- res[[3]]

save_plot2(p25, metabolite_biomarker_path, "randomforest_plot", width = 10, height = 8)

write_sheet2(metabolite_biomarker_wb, "randomforest_features", dat)
openxlsx::saveWorkbook(metabolite_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 29. loadingPCA.ms: 载荷矩阵挑选重要代谢物 ----------

res <- loadingPCA.ms(ps = ps.ms, Top = 20)
p <- res[[1]]
dat <- res[[2]]

save_plot2(p, metabolite_biomarker_path, "PCA_loading_plot", width = 10, height = 8)

write_sheet2(metabolite_biomarker_wb, "PCA_loading", dat)
openxlsx::saveWorkbook(metabolite_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 30. LDA.ms: LDA筛选特征代谢物 ----------

tablda <- LDA.ms(ps = ps.ms,
                 Top = 100,
                 p.lvl = 0.05,
                 lda.lvl = 1,
                 seed = 11,
                 adjust.p = FALSE)

p35 <- lefse_bar(taxtree = tablda[[2]])
dat <- tablda[[2]]

save_plot2(p35, metabolite_biomarker_path, "LDA_plot", width = 10, height = 8)

write_sheet2(metabolite_biomarker_wb, "LDA_features", dat)
openxlsx::saveWorkbook(metabolite_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 31. svm.ms: svm筛选特征微生物 ----------

res <- svm_ms(ps = pst %>% filter_OTU_ps(20), k = 5)
AUC <- res[[1]]
importance <- res[[2]]

write_sheet2(metabolite_biomarker_wb, "SVM_AUC", AUC)
write_sheet2(metabolite_biomarker_wb, "SVM_importance", importance)
openxlsx::saveWorkbook(metabolite_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 32. xgboost.ms: xgboost筛选特征微生物 ----------

library(xgboost)
library(Ckmeans.1d.dp)
library(mia)

res <- xgboost.ms(ps = pst, top = 20)
accuracy <- res[[1]]
importance <- res[[2]]
importance <- importance$importance

write_sheet2(metabolite_biomarker_wb, "XGBoost_accuracy", accuracy)
write_sheet2(metabolite_biomarker_wb, "XGBoost_importance", importance)
openxlsx::saveWorkbook(metabolite_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 33. glm.ms: glm筛选特征微生物 ----------

res <- glm.ms(ps = pst %>% filter_OTU_ps(50), k = 5)
AUC <- res[[1]]
importance <- res[[2]]

write_sheet2(metabolite_biomarker_wb, "GLM_AUC", AUC)
write_sheet2(metabolite_biomarker_wb, "GLM_importance", importance)
openxlsx::saveWorkbook(metabolite_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 34. lasso.ms: lasso筛选特征微生物 ----------

library(glmnet)

res <- lasso.ms(ps = pst, top = 20, seed = 1010, k = 5)
accuracy <- res[[1]]
importance <- res[[2]]

write_sheet2(metabolite_biomarker_wb, "Lasso_accuracy", accuracy)
write_sheet2(metabolite_biomarker_wb, "Lasso_importance", importance)
openxlsx::saveWorkbook(metabolite_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 35. decisiontree.ms: 决策树 ----------

library(rpart)

res <- decisiontree.ms(ps = ps.ms, top = 50, seed = 6358, k = 5)
accuracy <- res[[1]]
importance <- res[[2]]

write_sheet2(metabolite_biomarker_wb, "DecisionTree_accuracy", accuracy)
write_sheet2(metabolite_biomarker_wb, "DecisionTree_importance", importance)
openxlsx::saveWorkbook(metabolite_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 36. nnet.ms: 神经网络筛选特征微生物 ----------

res <- nnet.ms(ps = pst, top = 100, seed = 1010, k = 5)
accuracy <- res[[1]]
importance <- res[[2]]

write_sheet2(metabolite_biomarker_wb, "NeuralNet_accuracy", accuracy)
write_sheet2(metabolite_biomarker_wb, "NeuralNet_importance", importance)
openxlsx::saveWorkbook(metabolite_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 37. bagging.ms: Bootstrap Aggregating筛选特征微生物 ----------

library(ipred)

res <- bagging.ms(ps = pst, top = 20, seed = 1010, k = 5)
accuracy <- res[[1]]
importance <- res[[2]]

write_sheet2(metabolite_biomarker_wb, "Bagging_accuracy", accuracy)
write_sheet2(metabolite_biomarker_wb, "Bagging_importance", importance)
openxlsx::saveWorkbook(metabolite_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

## ===================== 7. 网络分析 =====================

metabolite_network_path <- file.path(metabolite_path, "06_network")
dir.create(metabolite_network_path, recursive = TRUE, showWarnings = FALSE)

network_xlsx_path <- file.path(metabolite_network_path, "network_results.xlsx")

if (file.exists(network_xlsx_path)) {
  metabolite_network_wb <- openxlsx::loadWorkbook(network_xlsx_path)
} else {
  metabolite_network_wb <- openxlsx::createWorkbook()
}

### ---------- 38. network.pip: 网络分析主函数 ----------

library(igraph)
detach(package:mia)
tax <- ps.ms %>% vegan_tax() %>% as.data.frame()
tax$ID <- NULL
tax_table(ps.ms) <- phyloseq::tax_table(as.matrix(tax))
head(tax)
tab.r <- network.pip(
  ps = ps.ms,
  N = 400,
  big = TRUE,
  select_layout = TRUE,
  layout_net = "model_maptree2",
  r.threshold = 0.6,
  p.threshold = 0.05,
  maxnode = 2,
  label = FALSE,
  lab = "elements",
  group = "Group",
  fill = "Super_class",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  step = 100,
  R = 10,
  ncpus = 1
)

dat <- tab.r[[2]]
cortab <- dat$net.cor.matrix$cortab

plot <- tab.r[[1]]
p0 <- plot[[1]]

save_plot2(p0, metabolite_network_path, "network_plot", width = 12, height = 10)

### ---------- 39. net_properties.4: 网络属性计算 ----------

i <- 1
cor <- cortab
id <- names(cor)

for (i in 1:length(id)) {
  igraph <- cor[[id[i]]] %>% make_igraph()
  dat <- net_properties.4(igraph, n.hub = FALSE)
  colnames(dat) <- id[i]

  if (i == 1) {
    dat2 <- dat
  } else{
    dat2 <- cbind(dat2, dat)
  }
}

write_sheet2(metabolite_network_wb, "network_properties", dat2)
openxlsx::saveWorkbook(metabolite_network_wb, network_xlsx_path, overwrite = TRUE)

### ---------- 40. netproperties.sample: 单个样本的网络属性 ----------

for (i in 1:length(id)) {
  ps.mst <- ps.ms %>% subset_samples.wt("Group", id[i]) %>% remove.zero()
  dat.f <- netproperties.sample(pst = ps.mst, cor = cor[[id[i]]])

  if (i == 1) {
    dat.f2 <- dat.f
  } else{
    dat.f2 <- rbind(dat.f2, dat.f)
  }
}

map <- sample_data(ps.ms)
map$ID <- row.names(map)
map <- map %>% as.tibble()
dat3 <- dat.f2 %>% rownames_to_column("ID") %>% inner_join(map, by = "ID")

write_sheet2(metabolite_network_wb, "sample_network_properties", dat3)
openxlsx::saveWorkbook(metabolite_network_wb, network_xlsx_path, overwrite = TRUE)

### ---------- 41. node_properties: 计算节点属性 ----------

for (i in 1:length(id)) {
  igraph <- cor[[id[i]]] %>% make_igraph()
  nodepro <- node_properties(igraph) %>% as.data.frame()
  nodepro$Group <- id[i]
  colnames(nodepro) <- paste0(colnames(nodepro), ".", id[i])
  nodepro <- nodepro %>%
    as.data.frame() %>%
    rownames_to_column("ASV.name")

  if (i == 1) {
    nodepro2 <- nodepro
  } else{
    nodepro2 <- nodepro2 %>% full_join(nodepro, by = "ASV.name")
  }
}

write_sheet2(metabolite_network_wb, "node_properties", nodepro2)
openxlsx::saveWorkbook(metabolite_network_wb, network_xlsx_path, overwrite = TRUE)

### ---------- 42. module.compare.net.pip: 网络显著性比较 ----------

dat <- module.compare.net.pip(
  ps = NULL,
  corg = cor,
  degree = TRUE,
  zipi = FALSE,
  r.threshold = 0.8,
  p.threshold = 0.05,
  method = "spearman",
  padj = FALSE,
  n = 3)
res <- dat[[1]]

write_sheet2(metabolite_network_wb, "network_comparison", res)
openxlsx::saveWorkbook(metabolite_network_wb, network_xlsx_path, overwrite = TRUE)

## ===================== 8. 通路富集分析 =====================

metabolite_pathway_path <- file.path(metabolite_path, "07_pathway")
dir.create(metabolite_pathway_path, recursive = TRUE, showWarnings = FALSE)

pathway_xlsx_path <- file.path(metabolite_pathway_path, "pathway_results.xlsx")

if (file.exists(pathway_xlsx_path)) {
  metabolite_pathway_wb <- openxlsx::loadWorkbook(pathway_xlsx_path)
} else {
  metabolite_pathway_wb <- openxlsx::createWorkbook()
}

if ("package:mia" %in% search()) {
  detach("package:mia", unload = TRUE)
}

### ---------- 43. pathway_enrich.ms: 通路富集 ----------
head(tax)

ps.ms3 <- ps.ms %>% tax_glom_wt("KEGG COMPOUND ID")

res2 <- pathway_enrich.ms(ps = ps.ms3, dif.method = "wilcox")

save_plot2(res2$plots$OE.WT.plot, metabolite_pathway_path, "pathway_enrich_OE_WT", width = 10, height = 8)

if(!is.null(res2$plots$KO.OE.plot)) {
  save_plot2(res2$plots$KO.OE.plot, metabolite_pathway_path, "pathway_enrich_KO_OE", width = 10, height = 8)
}

write_sheet2(metabolite_pathway_wb, "pathway_enrich_OE_WT", res2$plotdata$OE.WT)

if(!is.null(res2$plotdata$KO.OE)) {
  write_sheet2(metabolite_pathway_wb, "pathway_enrich_KO_OE", res2$plotdata$KO.OE)
}
openxlsx::saveWorkbook(metabolite_pathway_wb, pathway_xlsx_path, overwrite = TRUE)

### ---------- 44. reaction.show.ms: 反应展示 ----------

res3 <- reaction.show.ms(ps = ps.ms3, dif.method = "wilcox")

save_plot2(res3$plots$OE.WT.plot, metabolite_pathway_path, "reaction_show_OE_WT", width = 10, height = 8)

write_sheet2(metabolite_pathway_wb, "reaction_show_OE_WT", res3$plotdata$OE.WT)
openxlsx::saveWorkbook(metabolite_pathway_wb, pathway_xlsx_path, overwrite = TRUE)

### ---------- 45. bubble.enrich.ms: 气泡图展示富集分析结果 ----------

res4 <- buplotall.ms(ps = ps.ms3, dif.method = "wilcox")

save_plot2(res4$plots$OE.WT.plot, metabolite_pathway_path, "bubble_enrich_OE_WT", width = 10, height = 8)

write_sheet2(metabolite_pathway_wb, "bubble_enrich_OE_WT", res4$plotdata$OE.WT)
openxlsx::saveWorkbook(metabolite_pathway_wb, pathway_xlsx_path, overwrite = TRUE)




# 代谢组学特征挑选 -----

map = ps.ms %>%  sample_data()
map$Group
ps.ms2 = ps.ms %>% subset_taxa.wt("KEGG COMPOUND ID",NA,TRUE)
tax = ps.ms2 %>% vegan_tax() %>% as.data.frame()
head(tax)

tab_ms_selection = ps.ms2 %>% subset_samples.wt("Group",c("MM","MT")) %>% vegan_otu() %>%
  as.data.frame()

map = sample_data(ps.ms2 %>% subset_samples.wt("Group",c("MM","MT")))
map$Group = as.factor(map$Group)
map$ID = NULL
map
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
tax_ms = ps.ms2 %>% phyloseq::tax_table() %>% as.data.frame()
head(tax_ms)
tax_ms$metab_id = row.names(tax_ms)
dat_ms_final = dat_ms %>% left_join(tax_ms, by = c("Metabolite"="metab_id"))
we
head(dat_ms_final)
p = plot_physeq_metabolites(ps = ps.ms2 %>% subset_samples.wt("Group",c("MM","MT")),
                            dat_ms = dat_ms_final,top_n = 50)


p
ggsave("./result/metabolite/feature.pdf",width = 24,height = 10)




