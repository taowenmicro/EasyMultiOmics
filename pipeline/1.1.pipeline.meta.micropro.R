# 宏基因组物种分析#-----

rm(list = ls())

## 0. 基础设置 & 通用函数 ----------------------------------------------

# 加载包
# BiocManager::install("MicrobiotaProcess")
library(EasyMultiOmics)
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(ggrepel)
library(ggsci)
library(openxlsx)
library(ape)
library(picante)

# 结果输出路径（根据项目修改一次即可）
path <- "./result/metm/"
dir.create(path, showWarnings = FALSE, recursive = TRUE)

# 分组调色函数（可选，用于根据分组生成颜色）
get_group_cols <- function(groups,
                           palette = c("npg", "nejm", "lancet")) {
  palette <- match.arg(palette)
  
  pal_fun <- switch(
    palette,
    npg    = ggsci::pal_npg("nrc"),
    nejm   = ggsci::pal_nejm(),
    lancet = ggsci::pal_lancet()
  )
  
  cols <- pal_fun(length(groups))
  setNames(cols, groups)
}

# 通用 ggplot 保存函数
save_plot2 <- function(p,
                       out_dir,
                       prefix,
                       width       = NULL,   # 如果指定，就完全按这个宽高
                       height      = NULL,
                       dpi         = 300,
                       base_width  = 8,      # 未指定宽高时，从这个基础值自动放大
                       base_height = 6,
                       x_per_inch      = 6,  # 每 inch 放多少个 x 类别
                       facets_per_row  = 4   # 每行 facet 预期数量
) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 估算 x 轴类别数（如果拿不到就忽略）
  n_x <- NA_integer_
  if (!is.null(p$mapping$x)) {
    xvar <- rlang::as_name(p$mapping$x)
    if (!is.null(p$data) && xvar %in% names(p$data)) {
      n_x <- dplyr::n_distinct(p$data[[xvar]])
    }
  }
  
  # 估算 facet 数量
  n_facets <- 1L
  gb <- tryCatch(ggplot2::ggplot_build(p),
                 error = function(e) NULL)
  if (!is.null(gb) && !is.null(gb$layout$layout$PANEL)) {
    n_facets <- length(unique(gb$layout$layout$PANEL))
  }
  
  # 宽度：如果用户没指定，就根据 x 类别数加宽
  if (is.null(width)) {
    add_w <- if (!is.na(n_x)) max(0, n_x / x_per_inch) else 0
    width_final <- base_width + add_w
  } else {
    width_final <- width
  }
  
  # 高度：如果用户没指定，就根据 facet 行数加高
  if (is.null(height)) {
    facet_rows <- ceiling(n_facets / facets_per_row)
    add_h      <- max(0, (facet_rows - 1) * 2)  # 每多一行 facet 高 +2 inch
    height_final <- base_height + add_h
  } else {
    height_final <- height
  }
  
  # 实际保存 PNG + PDF
  ggplot2::ggsave(
    filename  = file.path(out_dir, paste0(prefix, ".png")),
    plot      = p,
    width     = width_final,
    height    = height_final,
    dpi       = dpi,
    units     = "in",
    limitsize = FALSE
  )
  
  ggplot2::ggsave(
    filename  = file.path(out_dir, paste0(prefix, ".pdf")),
    plot      = p,
    width     = width_final,
    height    = height_final,
    units     = "in"
  )
}

# circlize / base 图的保存函数
save_circlize_plot <- function(expr,
                               out_dir,
                               prefix,
                               width  = 10,
                               height = 10,
                               res    = 300) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 保存 PNG
  png(
    filename = file.path(out_dir, paste0(prefix, ".png")),
    width    = width,
    height   = height,
    units    = "in",
    res      = res
  )
  eval(substitute(expr), envir = parent.frame())
  dev.off()
  
  # 保存 PDF
  pdf(
    file   = file.path(out_dir, paste0(prefix, ".pdf")),
    width  = width,
    height = height
  )
  eval(substitute(expr), envir = parent.frame())
  dev.off()
}

# 统一写入 excel sheet 的函数
write_sheet2 <- function(wb, sheet_name, df, row_names = TRUE) {
  if (sheet_name %in% names(wb)) {
    openxlsx::removeWorksheet(wb, sheet_name)
  }
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, df, rowNames = row_names)
}

save_comp_wb <- function(wb, path) {
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  cat("已保存到：", normalizePath(path, winslash = "/"), "\n")
}

## 1. 宏基因组物种分析参数设定 -------------------------------------------

# 这里假定 ps.micro 已经在环境中，是一个宏基因组物种水平的 phyloseq 对象

# 有效基因数量
sample_sums(ps.micro)

# map 信息
map <- sample_data(ps.micro)
head(map)

phyloseq::tax_table(ps.micro) %>% head()

# 提取分组因子数量
gnum <- phyloseq::sample_data(ps.micro)$Group %>% unique() %>% length()
gnum

# 设定排序顺序：按照 ps.micro 对象中 map 文件的 Group 顺序
axis_order <- phyloseq::sample_data(ps.micro)$Group %>% unique()
axis_order

# 分组颜色（按你给的方案）
col.g <- c(
  "KO" = "#D55E00",
  "WT" = "#0072B2",
  "OE" = "#009E73"
)

# 主题设置（ggClusterNet / EasyMultiOmics）
package.amp()

res <- theme_my(ps.micro)
mytheme1 <- res[[1]]
mytheme2 <- res[[2]]
colset1  <- res[[3]]
colset2  <- res[[4]]
colset3  <- res[[5]]
colset4  <- res[[6]]

# alpha diversity -----
#1 alpha.metm: 6中alpha多样性计算#----
all.alpha = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )
#--alpha多样性指标运算
tab = alpha.metm(ps = ps.micro,group = "Group" )
head(tab)

data = cbind(data.frame(ID = 1:length(tab$Group),group = tab$Group),tab[all.alpha])
head(data)
data$ID = as.character(data$ID)
# data$Inv_Simpson[is.na(data$Inv_Simpson)]
# data$Inv_Simpson %>% tail(1000)

result = MuiKwWlx2(data = data,num = 3:6)

result1 = FacetMuiPlotresultBox(data = data,num = 3:6,
                                          result = result,
                                          sig_show ="abc",ncol = 4,width = 0.4 )


p1_1 = result1[[1]] +scale_fill_manual(values = col.g)

p1_1+
  ggplot2::scale_x_discrete(limits = axis_order) +
  # theme_cell()+
  theme_nature()+

  ggplot2::guides(fill = guide_legend(title = none))



res = FacetMuiPlotresultBar(data = data,num = c(3:6),
                                      result = result,sig_show ="abc",ncol = 4,
                            mult.y = 0.3
                            )
p1_2 = res[[1]]

p1_2 +
  ggplot2::scale_x_discrete(limits = axis_order) +
  theme_nature() +
  ggplot2::guides(fill = guide_legend(title = none)) +
  ggplot2::scale_fill_manual(values = col.g)

p = p1_2 +
  ggplot2::scale_x_discrete(limits = axis_order) +
  theme_cell() +
  ggplot2::guides(fill = guide_legend(title = none)) +
  ggplot2::scale_fill_manual(values = col.g)


res = FacetMuiPlotReBoxBar(data = data,num = c(3:6),result = result,sig_show ="abc",ncol = 4,
                           mult.y = 0.3,
                           lab.yloc = 1.1
                           )
p1_3 = res[[1]]
p1_3+
  ggplot2::scale_x_discrete(limits = axis_order) +
  theme_nature() +
  ggplot2::guides(fill = guide_legend(title = none)) +
  ggplot2::scale_fill_manual(values = col.g)

#基于输出数据使用ggplot出图
lab = result1[[2]] %>% distinct(group ,name,.keep_all = TRUE)
p1_0 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) +
  geom_violin(alpha=1, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
  labs(x="", y="")+
  facet_wrap(.~name,scales="free_y",ncol  = 4) +
  # theme_classic()+
  geom_text(data = lab,aes(x=group , y=y ,label=stat),alpha = 0.6)
p1_0+
  ggplot2::scale_x_discrete(limits = axis_order) +
  theme_nature() +
  ggplot2::guides(fill = guide_legend(title = none)) +
  ggplot2::scale_fill_manual(values = col.g)+
  ggplot2::scale_color_manual(values = col.g)

# 保存路径
alppath <- file.path(path, "alpha")
dir.create(alppath, showWarnings = FALSE)

# 保存Alpha多样性图片
save_plot2(
  p       = p1_1,
  out_dir = alppath,
  prefix  = "alpha_diversity_boxplot",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p1_2,
  out_dir = alppath,
  prefix  = "alpha_diversity_barplot",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p1_3,
  out_dir = alppath,
  prefix  = "alpha_diversity_reboxbar",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p1_0,
  out_dir = alppath,
  prefix  = "alpha_diversity_violinplot",
  width   = 10,
  height  = 8
)

#3 alpha_rare.metm:alpha多样性稀释曲线#---------
library(microbiome)
library(vegan)

rare <- mean(phyloseq::sample_sums(ps.micro))/10

result = alpha_rare.metm(ps = ps.micro ,
                          group = "Group", method = "Richness", start = 100, step = rare)

#--提供单个样本溪稀释曲线的绘制
p2_1 <- result[[1]] +scale_color_manual(values = col.g)+
  theme_nature()
p2_1
## 提供数据表格，方便输出
raretab <- result[[2]]
head(raretab)
#--按照分组展示稀释曲线
p2_2 <- result[[3]]+
  theme_nature()
p2_2
#--按照分组绘制标准差稀释曲线
p2_3 <- result[[4]]+
  theme_nature()
p2_3

# 保存路径
alppath <- file.path(path, "alpha")
dir.create(alppath, showWarnings = FALSE)

# 保存稀释曲线图片

save_plot2(
  p       = p2_1,
  out_dir = alppath,
  prefix  = "alpha_rarefaction_curve",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p2_2,
  out_dir = alppath,
  prefix  = "alpha_rarefaction_group_curve",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p2_3,
  out_dir = alppath,
  prefix  = "alpha_rarefaction_group_sd_curve",
  width   = 10,
  height  = 8
)

# beta diversity  -----
#4 ordinate.metm: 排序分析#----------
result = ordinate.metm(ps = ps.micro, group = "Group", dist = "bray",
                        method = "PCoA", Micromet = "anosim", pvalue.cutoff = 0.05)
p3_1 = result[[1]]
p3_1 +
  scale_fill_manual(values = col.g)+
  scale_color_manual(values = col.g,guide = F) +
  theme_nature()+
  theme(axis.title.y = element_text(angle = 90))

#带标签图形出图
p3_2 = result[[3]]
p3_2 +
  scale_fill_manual(values = col.g)+
  scale_color_manual(values = col.g,guide = F) +
  theme_nature()+
  theme(axis.title.y = element_text(angle = 90))


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

p3_3 +
  scale_fill_manual(values = col.g)+
  scale_color_manual(values = col.g,guide = F)+
  theme_cell()+
  theme(axis.title.y = element_text(angle = 90))

# 创建beta多样性主目录
betapath <- file.path(path, "beta")
dir.create(betapath, showWarnings = FALSE, recursive = TRUE)

#保存排序分析图片（用 save_plot2）
save_plot2(
  p       = p3_1,
  out_dir = betapath,
  prefix  = "pcoa_basic",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p3_2,
  out_dir = betapath,
  prefix  = "pcoa_labeled",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p3_3,
  out_dir = betapath,
  prefix  = "pcoa_refined",
  width   = 10,
  height  = 8
)

# 保存排序分析数据（用 write_sheet2 + save_comp_wb） 
beta_wb <- createWorkbook()

# 如果不需要行名可以改成 row_names = FALSE
write_sheet2(beta_wb, "plotdata", plotdata, row_names = TRUE)
write_sheet2(beta_wb, "cent",     cent,     row_names = TRUE)
write_sheet2(beta_wb, "segs",     segs,     row_names = TRUE)

# 实际保存文件（可按需修改文件名）
save_comp_wb(beta_wb, file.path(betapath, "beta_results.xlsx"))

#5 MicroTest.metm:群落水平差异检测-------
dat1 = MicroTest.metm(ps = ps.micro, Micromet = "adonis", dist = "bray")
dat1

# 保存数据
addWorksheet(beta_wb, "adonis_results")
writeData(beta_wb, "adonis_results", dat1)

#6 pairMicroTest.metm:两两分组群落水平差异检测-------
dat2 = pairMicroTest.metm(ps = ps.micro, Micromet = "MRPP", dist = "bray")
dat2

# 保存两两比较数据
write_sheet2(
  wb        = beta_wb,
  sheet_name = "pairwise_MRPP",
  df        = dat2,
  row_names = TRUE   # 若不需要行名可改为 FALSE
)

save_comp_wb(beta_wb, file.path(betapath, "beta_results.xlsx"))
#7 mantal.metm：群落差异检测普鲁士分析#------
result <- mantal.metm(ps = ps.micro,
                       method =  "spearman",
                       group = "Group",
                       ncol = gnum,
                       nrow = 1)
data <- result[[1]]
data
p3_7 <- result[[2]]

p3_7 +
  # scale_fill_manual(values = colset1)+
  # scale_color_manual(values = colset1,guide = F) +
  theme_cell

#  Mantel 分析图片（用 save_plot2）
save_plot2(
  p       = p3_7_final,
  out_dir = mantelpath,
  prefix  = "mantel_plot",
  width   = 10,
  height  = 8
)

# Mantel 分析数据（用 write_sheet2 + save_comp_wb）
# 往同一个 beta_wb 里加一个 sheet 存 Mantel 结果
write_sheet2(
  wb         = beta_wb,
  sheet_name = "mantel_results",
  df         = mantel_data,
  row_names  = TRUE   # 若不需要行名可改为 FALSE
)

# 最后统一保存工作簿
save_comp_wb(beta_wb, file.path(betapath, "beta_results.xlsx"))

#8 cluster.metm:样品聚类-----
res = cluster.metm (ps= ps.micro,
                     hcluter_method = "complete",
                     dist = "bray",
                     cuttree = 3,
                     row_cluster = TRUE,
                     col_cluster =  TRUE)

p4 = res[[1]] +
  scale_fill_manual(values = col.g)

p4
p4_1 = res[[2]]
p4_1
p4_2 = res[[3]]
p4_2
dat = res[[4]] # 聚类矩阵
head(dat)

# 创建聚类分析保存目录
# 创建聚类分析保存目录
clusterpath <- file.path(path, "cluster")
dir.create(clusterpath, showWarnings = FALSE, recursive = TRUE)

# 聚类分析图片（用 save_plot2）
save_plot2(
  p       = p4,
  out_dir = clusterpath,
  prefix  = "cluster_heatmap",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p4_1,
  out_dir = clusterpath,
  prefix  = "cluster_dendrogram",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p4_2,
  out_dir = clusterpath,
  prefix  = "cluster_scatter",
  width   = 10,
  height  = 8
)

# 聚类数据（用 write_sheet2 + save_comp_wb）
cluster_wb <- createWorkbook()

write_sheet2(
  wb         = cluster_wb,
  sheet_name = "cluster_matrix",
  df         = dat,
  row_names  = TRUE   # 和之前 rowNames = TRUE 一致
)

save_comp_wb(cluster_wb, file.path(clusterpath, "cluster_results.xlsx"))

# compositipn -----
# 创建组成分析保存目录
comppath <- file.path(path, "composition")
dir.create(comppath, showWarnings = FALSE, recursive = TRUE)

#9 Ven.Upset.metm: 用于展示共有、特有的OTU/ASV----
# 分组小于6时使用
res = Ven.Upset.metm(ps =  ps.micro,
                      group = "Group",
                      N = 0.5,
                      size = 3)
p10.1 = res[[1]]
p10.1
p10.2 = res[[2]]
p10.2
p10.3 = res[[4]]
p10.3
dat = res[[3]]
head(dat)

# 保存Venn/Upset图片
gsave_plot2(
  p       = p10.1,
  out_dir = comppath,
  prefix  = "venn_diagram",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p10.2,
  out_dir = comppath,
  prefix  = "upset_plot",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p10.3,
  out_dir = comppath,
  prefix  = "venn_additional",
  width   = 10,
  height  = 8
)

#Venn / Upset 数据（用 write_sheet2 + save_comp_wb)
comp_wb <- createWorkbook()

write_sheet2(
  wb         = comp_wb,
  sheet_name = "venn_upset_data",
  df         = dat,
  row_names  = TRUE   # 对应原来 rowNames = TRUE
)

save_comp_wb(comp_wb, file.path(comppath, "comp_results.xlsx"))
#10 ggVen.Upset.metm:用于展示共有、特有的OTU/ASV----
# 分组小于6时使用
res = ggVen.Upset.metm(ps = ps.micro,group = "Group")

# 保存ggVenn图片 (grid对象需要特殊处理)
png(file.path(comppath, "ggvenn_plot.png"), width = 10, height = 8, units = "in", res = 300)
grid::grid.draw(res[[1]])
dev.off()

pdf(file.path(comppath, "ggvenn_plot.pdf"), width = 10, height = 8)
grid::grid.draw(res[[1]])
dev.off()

dat = res[[2]]
dat
# 保存
write_sheet2(
  wb         = comp_wb,
  sheet_name = "ggvenn_data",
  df         = dat,
  row_names  = TRUE   # 对应原来 rowNames = TRUE
)

save_comp_wb(
  wb   = comp_wb,
  path = file.path(comppath, "comp_results.xlsx")
)
#11 VenSeper.micro: 详细展示每一组中OTU的物种组成 错-------
#---每个部分
library("ggpubr")
library(agricolae)
library(reshape2)

result = VenSuper.metm(ps = ps.micro,
                        group = "Group",
                        num = 6
)
# 提取韦恩图中全部部分的otu极其丰度做门类柱状图
p7_1 <- result[[1]]
p7_1+
  # scale_fill_manual(values = colset1)+
  # scale_color_manual(values = colset1,guide = F) +
  mytheme1
# 每部分的otu门类冲积图
p7_2 <- result[[3]]
p7_2+
  # scale_fill_manual(values = colset1)+
  # scale_color_manual(values = colset1,guide = F) +
  theme_bw()
# 保存
save_plot2(
  p       = p7_1,
  out_dir = comppath,
  prefix  = "vensuper_barplot",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p7_2,
  out_dir = comppath,
  prefix  = "vensuper_alluvial",
  width   = 10,
  height  = 8
)

#12 ggflower.metm:花瓣图展示共有特有微生物------
res <- ggflower.metm(ps = ps.micro ,
                       # rep = 1,
                       group = "ID",
                       start = 1, # 风车效果
                       m1 = 2, # 花瓣形状，方形到圆形到棱形，数值逐渐减少。
                       a = 0.3, # 花瓣胖瘦
                       b = 1, # 花瓣距离花心的距离
                       lab.leaf = 1, # 花瓣标签到圆心的距离
                       col.cir = "yellow",
                       N = 0.5 )

p13.1 = res[[1]]
p13.1

dat = res[[2]]
head(dat)

# 花瓣图图片（用 save_plot2）
save_plot2(
  p       = p13.1,
  out_dir = comppath,
  prefix  = "flower_plot",
  width   = 10,
  height  = 8
)

# 花瓣图数据（用 write_sheet2）
# 这里继续往前面已经创建的 comp_wb 里加一个 sheet
write_sheet2(
  wb         = comp_wb,
  sheet_name = "flower_data",
  df         = dat,
  row_names  = TRUE   # 对应原来 rowNames = TRUE
)

# 如果你还没有保存 comp_wb，可以在所有 sheet 都写完后统一保存一次，例如：
save_comp_wb(comp_wb, file.path(comppath, "comp_results.xlsx"))

#13 ven.network.metm:venn 网络----
result = ven.network.metm(
  ps = ps.micro ,#%>%  filter_OTU_ps(10000),
  N = 0.5,
  fill = "Phylum")
p14  = result[[1]]
p14

dat = result[[2]]
head(dat)

#Venn 网络图片（用 save_plot2）
save_plot2(
  p       = p14,
  out_dir = comppath,
  prefix  = "venn_network",
  width   = 10,
  height  = 8
)

# Venn 网络数据（用 write_sheet2，保留行名）
write_sheet2(
  wb         = comp_wb,
  sheet_name = "venn_network_data",
  df         = dat,
  row_names  = TRUE
)

# comp_wb 最后仍然用之前的 save_comp_wb 统一保存即可

 save_comp_wb(comp_wb, file.path(comppath, "comp_results.xlsx"))
#14 Micro_tern.metm: 三元图展示组成----
ps1 = ps.micro %>% filter_OTU_ps(500)
res = Micro_tern.metm(ps1)
p15 = res[[1]]
p15[[1]] +theme_classic()

dat =  res[[2]]
head(dat)

# 三元图图片（用 save_plot2）
save_plot2(
  p       = p15[[1]] + theme_classic(),
  out_dir = comppath,
  prefix  = "ternary_plot",
  width   = 10,
  height  = 8
)

# 三元图数据（用 write_sheet2，保留行名） 
write_sheet2(
  wb         = comp_wb,
  sheet_name = "ternary_data",
  df         = dat,
  row_names  = TRUE
)

# comp_wb 依然在最后统一用 save_comp_wb 保存即可，
 save_comp_wb(comp_wb, file.path(comppath, "comp_results.xlsx"))
#15 barMainplot.metm: 堆积柱状图展示组成----
library(ggalluvial)
# pst = ps.micro %>% subset_taxa.wt("Species","Unassigned",TRUE)
pst = ps.micro %>% subset_taxa.wt("Genus","Unassigned",TRUE)

result = barMainplot.metm(ps = pst,
                           j = "Genus",
                           # axis_ord = axis_order,
                           label = FALSE,
                           sd = FALSE,
                           Top =5)

 p4_1 <- result[[1]]
# scale_fill_manual(values = colset2) +
# scale_x_discrete(limits = axis_order) +
# mytheme1
p4_1+
  scale_fill_manual(values = colset2) +
  scale_x_discrete(limits = axis_order) +
  theme_nature()+
  theme(axis.title.y = element_text(angle = 90))


p4_2  <- result[[3]]
# # scale_fill_brewer(palette = "Paired") +
# scale_fill_manual(values = colset2) +
# scale_x_discrete(limits = axis_order) +
# mytheme1
p4_2+
  scale_fill_manual(values = colset2) +
  scale_x_discrete(limits = axis_order) +
  theme_nature()+
  theme(axis.title.y = element_text(angle = 90))

databar <- result[[2]] %>% group_by(Group,aa) %>%
  dplyr::summarise(sum(Abundance)) %>% as.data.frame()
colnames(databar) = c("Group","genus","Abundance(%)")
head(databar)

# 堆积柱状图图片（用 save_plot2）
  save_plot2(
    p       = p4_1,
    out_dir = comppath,
    prefix  = "barplot_samples",
    width   = 10,
    height  = 8
  )

save_plot2(
  p       = p4_2,
  out_dir = comppath,
  prefix  = "barplot_groups",
  width   = 10,
  height  = 8
)

# 堆积柱状图数据（用 write_sheet2）
# 原来：rowNames = TRUE
write_sheet2(
  wb         = comp_wb,
  sheet_name = "barplot_raw_data",
  df         = result[[2]],
  row_names  = TRUE
)

# 汇总数据通常不需要行名
write_sheet2(
  wb         = comp_wb,
  sheet_name = "barplot_summary_data",
  df         = databar,
  row_names  = FALSE
)

# 工作簿 comp_wb 仍然在最后统一用 save_comp_wb 保存即可，例如：
save_comp_wb(comp_wb, file.path(comppath, "comp_results.xlsx"))

#16 cluMicro.bar.metm: 聚类堆积柱状图展示组成-----
result <-  cluMicro.bar.metm (dist = "bray",
                               ps = ps.micro,
                               j = "Genus",
                               Top = 7, # 提取丰度前十的物种注释
                               tran = TRUE, # 转化为相对丰度值
                               hcluter_method = "complete",
                               Group = "Group",
                               cuttree = length(unique(phyloseq::sample_data(ps.micro)$Group)))

p5_2 <- result[[2]]
p5_2

p5_3 <- result[[3]]
p5_3
p5_4 <- result[[4]]
p5_4
clubardata <- result[[5]]

# 聚类堆积柱状图图片（用 save_plot2）
save_plot2(
  p       = p5_2,
  out_dir = comppath,
  prefix  = "cluster_barplot_1",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p5_3,
  out_dir = comppath,
  prefix  = "cluster_barplot_2",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p5_4,
  out_dir = comppath,
  prefix  = "cluster_barplot_3",
  width   = 10,
  height  = 8
)

# 聚类堆积柱状图数据（用 write_sheet2）
write_sheet2(
  wb         = comp_wb,
  sheet_name = "cluster_barplot_data",
  df         = clubardata,
  row_names  = TRUE   # 对应原来 rowNames = TRUE
)

# comp_wb 依旧在所有 sheet 写完后，用 save_comp_wb 统一保存即可：
save_comp_wb(comp_wb, file.path(comppath, "comp_results.xlsx"))
#17 cir_barplot.metf:环状堆积柱状图 -----
library(ggtree) # j = "Phylum"
res = cir_barplot.metm(
  ps = ps.micro,
  Top = 7,
  dist = "bray",
  cuttree = 3,
  hcluter_method = "complete")

p17 = res[[1]]
p17
dat= res[[2]]
head(dat)

# 环状堆积柱状图图片（用 save_plot2）
save_plot2(
  p       = p17,
  out_dir = comppath,
  prefix  = "circular_barplot",
  width   = 10,
  height  = 8
)

#  环状堆积柱状图数据（用 write_sheet2，保留行名） 
write_sheet2(
  wb         = comp_wb,
  sheet_name = "circular_barplot_data",
  df         = dat,
  row_names  = TRUE
)

# comp_wb 仍在最后统一用 save_comp_wb 保存，例如：
save_comp_wb(comp_wb, file.path(comppath, "comp_results.xlsx"))
#18 cir_plot.metm:和弦图展示物种组成-----
res = cir_plot.metm(ps  = ps.micro,Top = 12,rank = 6)
p <- recordPlot()

# 和弦图图片（用 save_circlize_plot） 
  save_circlize_plot(
    expr   = replayPlot(p),   # 或 expr = { replayPlot(p) }
    out_dir = comppath,
    prefix  = "cir_plot",
    width   = 10,
    height  = 8,
    res     = 300
  )

# 和弦图数据（用 write_sheet2，保留行名）
write_sheet2(
  wb         = comp_wb,
  sheet_name = "cir_plot_data",
  df         = res[[1]],
  row_names  = TRUE
)

# comp_wb 仍在所有 sheet 写完后统一保存，例如：
save_comp_wb(comp_wb, file.path(comppath, "comp_results.xlsx"))
#19 Microheatmap.metm: 热图展示物种相对丰度差异-----
heatnum = 30
# map = phyloseq::sample_data(ps)
# map$ID = row.names(map)
# phyloseq::sample_data(ps) = map

ps_tem = ps.micro %>%
  ggClusterNet::scale_micro(method = "TMM") %>%
  ggClusterNet::tax_glom_wt(ranks = "Genus")
id <- ps.micro %>%
  ggClusterNet::scale_micro(method = "TMM") %>%
  ggClusterNet::tax_glom_wt(ranks = "Genus") %>%
  ggClusterNet::filter_OTU_ps(100) %>%
  ggClusterNet::vegan_otu() %>%
  t() %>% as.data.frame() %>%rowCV %>%
  sort(decreasing = TRUE) %>%
  head(heatnum) %>%
  names()
#? Microheatmap.metm
result <- Microheatmap.metm(ps_rela = ps_tem,id = id ,col_cluster = FALSE,
                            row_cluster = FALSE,
                            col1 = (ggsci::pal_gsea(alpha = 1))(12))

p24.1 <- result[[1]]

p24.1
# p1 +  scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(60))
p24.2 <- result[[2]]
p24.2
dat = result[[3]]
head(dat)

#  微生物热图图片（用 save_plot2） 
save_plot2(
  p       = p24.1,
  out_dir = comppath,
  prefix  = "microbiome_heatmap_1",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p24.2,
  out_dir = comppath,
  prefix  = "microbiome_heatmap_2",
  width   = 10,
  height  = 8
)

# 微生物热图数据（用 write_sheet2）
write_sheet2(
  wb         = comp_wb,
  sheet_name = "heatmap_data",
  df         = dat,
  row_names  = TRUE      # 保留行名
)

write_sheet2(
  wb         = comp_wb,
  sheet_name = "selected_otus",
  df         = data.frame(OTU_ID = id),
  row_names  = FALSE     # 原来未指定 rowNames，默认不写行名
)

# 统一保存 comp_wb
save_comp_wb(comp_wb, file.path(comppath, "comp_results.xlsx"))
# difference analysis -----
# 创建差异分析主目录
diffpath <- file.path(path, "differential")
dir.create(diffpath, showWarnings = FALSE, recursive = TRUE)
#20 EdgerSuper.metm:EdgeR计算差异微生物----
res = EdgerSuper.metm(ps = ps.micro %>% ggClusterNet::filter_OTU_ps(500),group  = "Group",artGroup = NULL, j = "Species")

p25.1 =  res[[1]][1]
p25.1
p25.2 =  res[[1]][2]
p25.2
p25.3 =  res[[1]][3]
p25.3
dat =  res[[2]]
head(dat)

#  EdgeR 图片（用 save_plot2）
save_plot2(
  p       = p25.1,
  out_dir = diffpath,
  prefix  = "edger_volcano_1",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p25.2,
  out_dir = diffpath,
  prefix  = "edger_volcano_2",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p25.3,
  out_dir = diffpath,
  prefix  = "edger_volcano_3",
  width   = 10,
  height  = 8
)

# EdgeR 数据（用 write_sheet2 + save_comp_wb）
diff_wb <- createWorkbook()

write_sheet2(
  wb         = diff_wb,
  sheet_name = "edger_results",
  df         = dat,
  row_names  = TRUE      # 对应原来 rowNames = TRUE
)

save_comp_wb(
  wb   = diff_wb,
  path = file.path(diffpath, "diff_results.xlsx")
)
#21 EdgerSuper2.metm:EdgeR计算差异微生物-----
res =  EdgerSuper2.metm (ps = ps.micro,group  = "Group",artGroup =NULL, j = "Species")
head(res)

# 保存EdgeR2数据
write_sheet2(
  wb         = diff_wb,
  sheet_name = "edger2_results",
  df         = res,
  row_names  = TRUE
)

# 之后统一保存 diff_wb 即可，例如：
save_comp_wb(diff_wb, file.path(diffpath, "diff_results.xlsx"))
#22 DESep2Super.metm:DESep2计算差异微生物-----
res = DESep2Super.metm(ps = ps.micro %>% ggClusterNet::filter_OTU_ps(500),
                        group  = "Group",
                        artGroup =NULL,
                        j = "Species"
                        # path = diffpath.1
)

p26.1 =  res[[1]][1]
p26.1
p26.2 =  res[[1]][2]
p26.2
p26.3 =  res[[1]][3]
p26.3
dat =  res[[2]]
dat

# DESeq2 图片（用 save_plot2）
save_plot2(
  p       = p26.1[[1]],
  out_dir = diffpath,
  prefix  = "deseq2_volcano_1",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p26.2[[1]],
  out_dir = diffpath,
  prefix  = "deseq2_volcano_2",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p26.3[[1]],
  out_dir = diffpath,
  prefix  = "deseq2_volcano_3",
  width   = 10,
  height  = 8
)

# DESeq2 数据（用 write_sheet2）
write_sheet2(
  wb         = diff_wb,
  sheet_name = "deseq2_results",
  df         = dat,
  row_names  = TRUE    # 对应原来 rowNames = TRUE
)

# diff_wb 最后统一保存，例如：
save_comp_wb(diff_wb, file.path(diffpath, "diff_results.xlsx"))
#23 edge_Manhattan.metm: 曼哈顿图展示差异微生物------
res = edge_Manhattan.metm(
  ps = ps.micro%>% ggClusterNet::filter_OTU_ps(500),
  pvalue = 0.05,
  lfc = 0
)
p27.1= res[[1]]
p27.1
p27.2= res[[2]]
p27.2
p27.3= res[[3]]
p27.3
# 保存
save_plot2(
  p       = p27.1,
  out_dir = diffpath,
  prefix  = "manhattan_plot_1",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p27.2,
  out_dir = diffpath,
  prefix  = "manhattan_plot_2",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p27.3,
  out_dir = diffpath,
  prefix  = "manhattan_plot_3",
  width   = 10,
  height  = 8
)
#24 stemp_diff.metm: stamp展示差异微生物----
# map = phyloseq::sample_data(ps)
# map$ID = row.names(map)
#
# # map$Group = as.factor(map$Group)
# sample_data(ps) = map
allgroup <- combn(unique(map$Group),2)
plot_list <- list()
for (i in 1:dim(allgroup)[2]) {
  ps_sub <- phyloseq::subset_samples(ps.micro,Group %in% allgroup[,i]);ps_sub
  p <- stemp_diff.metm(ps = ps_sub,Top = 20,ranks = 6)
  plot_list[[i]] <- p

}

p28.1 = plot_list[[1]]
p28.1
p28.2 = plot_list[[2]]
p28.2
p28.3 = plot_list[[3]]
p28.3

#  STAMP 图图片（用 save_plot2）
save_plot2(
  p       = p28.1,
  out_dir = diffpath,
  prefix  = "stamp_plot_1",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p28.2,
  out_dir = diffpath,
  prefix  = "stamp_plot_2",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p28.3,
  out_dir = diffpath,
  prefix  = "stamp_plot_3",
  width   = 10,
  height  = 8
)
#25 Mui.Group.volcano.metm: 聚类火山图------
res =  EdgerSuper2.metm (ps = ps.micro,group  = "Group",artGroup =NULL, j = "OTU")
res2 = Mui.Group.volcano.metm(res = res)

p29.1 = res2[[1]]
p29.1
p29.2 = res2[[2]]
p29.2
dat = res2[[3]]
dat

# 聚类火山图图片（用 save_plot2） 
save_plot2(
  p       = p29.1,
  out_dir = diffpath,
  prefix  = "cluster_volcano_1",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p29.2,
  out_dir = diffpath,
  prefix  = "cluster_volcano_2",
  width   = 10,
  height  = 8
)

# 聚类火山图数据（用 write_sheet2 + save_comp_wb）
write_sheet2(
  wb         = diff_wb,
  sheet_name = "cluster_volcano_results",
  df         = dat,
  row_names  = TRUE
)

save_comp_wb(
  wb   = diff_wb,
  path = file.path(diffpath, "diff_results.xlsx")
)
# biomarker identification -----
# 创建生物标记物分析主目录
biomarkerpath <- file.path(path, "biomarker")
dir.create(biomarkerpath, showWarnings = FALSE, recursive = TRUE)
#26 rfcv.metm :交叉验证结果-------
library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
library(e1071)
result = rfcv.metm(ps = ps.micro %>% filter_OTU_ps(100),
                   group  = "Group",optimal = 20,nrfcvnum = 6)

prfcv = result[[1]]

prfcv+theme_classic()

# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable

# 随机森林交叉验证图片（用 save_plot2）
save_plot2(
  p       = prfcv,
  out_dir = biomarkerpath,
  prefix  = "rfcv_plot",
  width   = 10,
  height  = 8
)

#  随机森林交叉验证数据（用 write_sheet2 + save_comp_wb）
biomarker_wb <- createWorkbook()

write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "rfcv_plot_data",
  df         = result[[2]],
  row_names  = TRUE
)

write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "rfcv_summary_table",
  df         = rfcvtable,
  row_names  = TRUE
)

# 统一保存
save_comp_wb(
  wb   = biomarker_wb,
  path = file.path(biomarkerpath, "biomarker_results.xlsx")
)
#27 Roc.metm:ROC 曲线绘制----
id = sample_data(ps.micro)$Group %>% unique()
aaa = combn(id,2)
i= 1
group = c(aaa[1,i],aaa[2,i])
b= data.frame(group)

ps = ps.micro %>% subset_taxa.wt("Family","Unassigned",TRUE)
ps = ps %>% subset_taxa.wt("Order","Unassigned",TRUE)
ps = ps %>% subset_taxa.wt("Genus","Unassigned",TRUE)
ps = ps %>% subset_taxa.wt("Phylum","Unassigned",TRUE)
ps = ps %>% subset_taxa.wt("class","Unassigned",TRUE)

pst = ps %>% subset_samples.wt("Group",group) %>%
  filter_taxa(function(x) sum(x ) > 10, TRUE)
res = Roc.metm( ps = pst %>% filter_OTU_ps(1000),group  = "Group",repnum = 5)
p33.1 =  res[[1]]
p33.1+theme_classic()
AUC =  res[[2]]
AUC
dat =  res[[3]]
dat

#  ROC 分析图片（用 save_plot2）
save_plot2(
  p       = p33.1,
  out_dir = biomarkerpath,
  prefix  = "roc_curve",
  width   = 10,
  height  = 8
)

#  ROC 分析数据（用 write_sheet2） 
write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "roc_results",
  df         = dat,
  row_names  = TRUE
)

write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "auc_summary",
  df         = data.frame(AUC = AUC),
  row_names  = TRUE
)

# 所有 biomarker 结果写完后再统一保存，例如：
 save_comp_wb(biomarker_wb, file.path(biomarkerpath, "biomarker_results.xlsx"))
#28 loadingPCA.metm: 载荷矩阵筛选特征微生物------
res = loadingPCA.metm(ps = ps.micro,Top = 20)
p34.1 = res[[1]]
p34.1+theme_classic()
dat = res[[2]]
dat

#  PCA 载荷分析图片（用 save_plot2） 
save_plot2(
  p       = p34.1,
  out_dir = biomarkerpath,
  prefix  = "pca_loading",
  width   = 10,
  height  = 8
)

# PCA 载荷分析数据（用 write_sheet2）
write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "pca_loading_data",
  df         = dat,
  row_names  = TRUE
)

# biomarker_wb 依旧在所有 sheet 写完后统一保存，例如：
save_comp_wb(biomarker_wb, file.path(biomarkerpath, "biomarker_results.xlsx"))
#29 LDA.metm: LDA筛选特征微生物-----
tablda = LDA.metm(ps = ps.micro,
                   Top = 10,
                   p.lvl = 0.05,
                   lda.lvl = 4,
                   seed = 11,
                   adjust.p = F)

p35 <- lefse_bar(taxtree = tablda[[2]])
p35+theme_classic()
dat = tablda[[2]]
dat

# LDA 分析图片（用 save_plot2）
save_plot2(
  p       = p35,
  out_dir = biomarkerpath,
  prefix  = "lda_barplot",
  width   = 10,
  height  = 8
)

# LDA 分析数据（用 write_sheet2） 
write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "lda_results",
  df         = dat,
  row_names  = TRUE
)

write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "lda_parameters",
  df         = data.frame(
    Parameter = c("Top", "p.lvl", "lda.lvl", "seed", "adjust.p"),
    Value     = c("10", "0.05", "4", "11", "FALSE")
  ),
  row_names  = FALSE
)

# biomarker_wb 在所有 sheet 写完后统一保存，例如：
save_comp_wb(biomarker_wb, file.path(biomarkerpath, "biomarker_results.xlsx"))
#30 svm.metm:svm筛选特征微生物 ----
res <- svm.metm(ps = pst %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#  SVM 分析数据（用 write_sheet2）
write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "svm_auc",
  df         = data.frame(AUC = AUC),
  row_names  = FALSE   # 原来没写 rowNames，默认不需要
)

write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "svm_feature_importance",
  df         = importance,
  row_names  = TRUE    # 保留行名
)

# 所有 biomarker 结果写完后统一保存，例如：
save_comp_wb(biomarker_wb, file.path(biomarkerpath, "biomarker_results.xlsx"))
#31 glm.metm :glm筛选特征微生物----
res <- glm.metm (ps = pst %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#  GLM 分析数据（用 write_sheet2） 
write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "glm_auc",
  df         = data.frame(AUC = AUC),
  row_names  = FALSE   # 原来未写 rowNames，默认不需要
)

write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "glm_feature_importance",
  df         = importance,
  row_names  = TRUE    # 保留行名
)

# 所有 biomarker 结果写完后统一保存，例如：
save_comp_wb(biomarker_wb, file.path(biomarkerpath, "biomarker_results.xlsx"))
#32 xgboost.metm: xgboost筛选特征微生物----
library(xgboost)
library(Ckmeans.1d.dp)
library(mia)
res = xgboost.metm(ps =pst, top = 20  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance
importance_data <- as.data.frame(importance_data)

# XGBoost 分析数据（用 write_sheet2）
write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "xgboost_accuracy",
  df         = data.frame(Accuracy = accuracy),
  row_names  = FALSE   # 原来未指定 rowNames，默认不需要
)

write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "xgboost_feature_importance",
  df         = importance_data,
  row_names  = TRUE    # 保留行名
)

# 所有 biomarker 结果写完后统一保存，例如：
save_comp_wb(biomarker_wb, file.path(biomarkerpath, "biomarker_results.xlsx"))
#33 lasso.metm: lasso筛选特征微生物----
library(glmnet)
res =lasso.metm(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# Lasso 分析数据（用 write_sheet2）
write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "lasso_accuracy",
  df         = data.frame(Accuracy = accuracy),
  row_names  = FALSE   # 原来未指定 rowNames，默认不需要
)

write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "lasso_feature_importance",
  df         = importance,
  row_names  = TRUE    # 保留行名
)

# 所有 biomarker 结果写完后统一保存，例如：
save_comp_wb(biomarker_wb, file.path(biomarkerpath, "biomarker_results.xlsx"))
#34 decisiontree.micro: 错----
library(rpart)
res =decisiontree.metm(ps=pst, top = 50, seed = 6358, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#  决策树分析数据（用 write_sheet2）
write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "tree_accuracy",
  df         = data.frame(Accuracy = accuracy),
  row_names  = FALSE   # 原来未指定 rowNames，默认不需要
)

write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "tree_feature_importance",
  df         = importance,
  row_names  = TRUE    # 保留行名
)

# 所有 biomarker 结果写完后统一保存，例如：
save_comp_wb(biomarker_wb, file.path(biomarkerpath, "biomarker_results.xlsx"))
#35 naivebayes.metm: bayes筛选特征微生物----
res = naivebayes.metm(ps=pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#  朴素贝叶斯分析数据（用 write_sheet2）
write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "naivebayes_accuracy",
  df         = data.frame(Accuracy = accuracy),
  row_names  = FALSE   # 原来未指定 rowNames，默认不需要
)

write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "naivebayes_feature_importance",
  df         = importance,
  row_names  = TRUE    # 保留行名
)

# 所有 biomarker 结果写完后统一保存，例如：
 save_comp_wb(biomarker_wb, file.path(biomarkerpath, "biomarker_results.xlsx"))
#36 randomforest.metm: 随机森林筛选特征微生物----
res = randomforest.metm( ps = pst%>% filter_OTU_ps(20),group  = "Group", optimal = 20)
p42.1 = res[[1]]
p42.1+theme_classic()
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
p42.4 =res[[4]]
p42.4

# 随机森林图片（用 save_plot2） 
save_plot2(
  p       = p42.1,
  out_dir = biomarkerpath,
  prefix  = "rf_importance",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p42.2,
  out_dir = biomarkerpath,
  prefix  = "rf_error",
  width   = 10,
  height  = 8
)

save_plot2(
  p       = p42.4,
  out_dir = biomarkerpath,
  prefix  = "rf_additional",
  width   = 10,
  height  = 8
)

# 随机森林数据（用 write_sheet2）
write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "rf_feature_importance",
  df         = dat,
  row_names  = TRUE
)

# biomarker_wb 仍在所有 sheet 写完后统一保存，例如：
 save_comp_wb(biomarker_wb, file.path(biomarkerpath, "biomarker_results.xlsx"))
#37 bagging.metm: Bootstrap Aggregating筛选特征微生物 ------
library(ipred)
res =bagging.metm(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Bagging分析数据
addWorksheet(biomarker_wb, "bagging_accuracy")
addWorksheet(biomarker_wb, "bagging_feature_importance")
writeData(biomarker_wb, "bagging_accuracy", data.frame(Accuracy = accuracy))
writeData(biomarker_wb, "bagging_feature_importance", importance, rowNames = TRUE)

#38 nnet.metm : 神经网络筛选特征微生物 ------
res =nnet.metm (ps=pst, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存神经网络分析数据

write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "nnet_accuracy",
  df         = data.frame(Accuracy = accuracy),
  row_names  = FALSE   # 原来未指定 rowNames，默认不需要
)

write_sheet2(
  wb         = biomarker_wb,
  sheet_name = "nnet_feature_importance",
  df         = importance,
  row_names  = TRUE    # 保留行名
)

save_comp_wb(
  wb   = biomarker_wb,
  path = file.path(biomarkerpath, "biomarker_results.xlsx")
)
#network analysis -----
# 创建网络分析主目录
networkpath <- file.path(path, "network")
dir.create(networkpath, showWarnings = FALSE, recursive = TRUE)
#39 network.pip:网络分析主函数--------
library(ggClusterNet)
library(igraph)
detach("package:mia", unload = TRUE)
tab.r = network.pip(
  ps = ps.micro,
  N = 200,
  # ra = 0.05,
  big = TRUE,
  select_layout = FALSE,
  layout_net = "model_maptree2",
  r.threshold = 0.6,
  p.threshold = 0.05,
  maxnode = 2,
  # method = "sparcc",
  label = FALSE,
  lab = "elements",
  group = "Group",
  fill = "Phylum",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 100,
  R=10,
  ncpus = 1
)

dat = tab.r[[2]]

cortab = dat$net.cor.matrix$cortab

#-提取全部图片的存储对象
plot = tab.r[[1]]

# 提取网络图可视化结果
p0 = plot[[1]]
p0


# 保存相关性矩阵信息
if (!is.null(cortab)) {
  for (i in 1:length(cortab)) {
    group_name <- names(cortab)[i]
    addWorksheet(main_network_wb, paste0("correlation_", group_name))

    # 转换相关性矩阵为数据框
    if (is.matrix(cortab[[i]])) {
      cor_df <- as.data.frame(cortab[[i]])
      writeData(main_network_wb, paste0("correlation_", group_name), cor_df, rowNames = TRUE)
    }
  }
}

# 保存网络基本信息
network_info <- data.frame(
  Parameter = c("N_features", "r_threshold", "p_threshold", "maxnode", "layout", "cluster_method"),
  Value = c("200", "0.6", "0.05", "2", "model_maptree2", "cluster_fast_greedy")
)
# 主网络分析图片（用 save_plot2） 
save_plot2(
  p       = p0,
  out_dir = networkpath,
  prefix  = "main_network_plot",
  width   = 12,
  height  = 10
)

#  主网络分析数据（用 write_sheet2） 
network_wb <- createWorkbook()

write_sheet2(
  wb         = network_wb,
  sheet_name = "network_parameters",
  df         = network_info,
  row_names  = FALSE   # 原来未指定 rowNames，默认不需要
)

# 后续如果还有其它网络相关 sheet，可继续 write_sheet2(...)
# 全部写完后统一保存，例如：
save_comp_wb(network_wb, file.path(networkpath, "network_results.xlsx"))

#40 net_properties.4:网络属性计算#---------
i = 1
cor = cortab
id = names(cor)
for (i in 1:length(id)) {
  igraph= cor[[id[i]]] %>% make_igraph()
  dat = net_properties.4(igraph,n.hub = F)
  head(dat,n = 16)
  colnames(dat) = id[i]

  if (i == 1) {
    dat2 = dat
  } else{
    dat2 = cbind(dat2,dat)
  }
}
head(dat2)

# 保存网络属性数据
write_sheet2(
  wb         = network_wb,
  sheet_name = "network_properties_summary",
  df         = dat2,
  row_names  = TRUE    # 保留行名
)

# 后面所有网络相关 sheet 写完后统一保存，例如：
save_comp_wb(network_wb, file.path(networkpath, "network_results.xlsx"))
#41 netproperties.sample:单个样本的网络属性#-------
for (i in 1:length(id)) {
  pst = ps.micro %>% subset_samples.wt("Group",id[i]) %>% remove.zero()
  dat.f = netproperties.sample(pst = pst,cor = cor[[id[i]]])
  # head(dat.f)
  if (i == 1) {
    dat.f2 = dat.f
  } else{
    dat.f2 = rbind(dat.f2,dat.f)
  }
}

map = map %>% as.tibble()
dat3 = dat.f2 %>% rownames_to_column("ID") %>% inner_join(map,by = "ID")

head(dat3)

# 保存样本网络属性数据
write_sheet2(
  wb         = network_wb,
  sheet_name = "sample_network_properties",
  df         = dat.f2,
  row_names  = TRUE    # 原来 rowNames = TRUE
)

write_sheet2(
  wb         = network_wb,
  sheet_name = "sample_metadata_combined",
  df         = dat3,
  row_names  = FALSE   # 原来未指定 rowNames，默认不需要
)

# 所有网络相关 sheet 写完后统一保存，例如：
save_comp_wb(network_wb, file.path(networkpath, "network_results.xlsx"))
#42 node_properties:计算节点属性#---------
for (i in 1:length(id)) {
  igraph= cor[[id[i]]] %>% make_igraph()
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

# 保存节点属性数据
write_sheet2(
  wb         = network_wb,
  sheet_name = "node_properties_combined",
  df         = nodepro2,
  row_names  = FALSE   # 与原来 rowNames = FALSE 保持一致
)

# 所有网络相关 sheet 写完后统一保存，例如：
save_comp_wb(network_wb, file.path(networkpath, "network_results.xlsx"))

#43 module.compare.net.pip:网络显著性比较#-----
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
head(res)

# 保存网络比较数据
addWorksheet(network_wb, "network_comparison_results")
writeData(network_wb, "network_comparison_results", res)
saveWorkbook(network_wb, file.path(networkpath, "network_results.xlsx"), overwrite = TRUE)


