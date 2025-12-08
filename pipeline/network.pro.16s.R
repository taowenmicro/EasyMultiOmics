
# 两组细菌网络比对流程

rm(list=ls())
library(EasyMultiOmics)
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(EasyMultiOmics)
library(igraph)

psg = ps.16s %>%
  subset_samples.wt("Group",c("WT","OE")) %>%
  tax_glom(taxrank = "Genus", NArm = TRUE) %>%
  scale_micro()



map= sample_data(psg)
head(map)


# 提取分组因子数量
gnum = phyloseq::sample_data(ps.16s)$Group %>% unique() %>% length()
gnum
#--设定排序顺序1：按照ps.16s对象中map文件顺序进行
axis_order =  phyloseq::sample_data(psg)$Group %>%unique();axis_order

col.g <- get_group_cols(axis_order)
scales::show_col(col.g)

# 创建扩增子微生物组分析目录
amplicon_network_path =  "../result/network.bac/"
# fs::dir_create(amplicon_network_path)
dir.create(amplicon_network_path, showWarnings = FALSE, recursive = TRUE)


#-主题--颜色等
package.amp()
res = theme_my(psg)
mytheme1 = res[[1]]
mytheme2 = res[[2]];
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]



#0---network.pip:网络分析主函数--------
tax = psg %>% vegan_tax() %>% as.data.frame()
head(tax)

tab.r = network.pip(
  ps = psg,
  N = 0,
  # ra = 0.05,
  big = TRUE,
  select_layout = FALSE,
  layout_net = "model_maptree2",
  r.threshold = 0.6,
  p.threshold = 0.05,
  maxnode = 2,
  method = "spearman",
  label =FALSE,
  lab = "elements",
  group = "Group",
  fill = "Phylum",
  size = "igraph.degree",
  zipi = F,
  ram.net = F,
  clu_method = "cluster_fast_greedy",
  step = 100,
  R=10,
  ncpus = 1
)


dat = tab.r[[2]]
cortab = dat$net.cor.matrix$cortab



# 1--保存相关矩阵 csv#-------------
for (nm in names(cortab)) {
  mat <- cortab[[nm]]
  out_file <- file.path(amplicon_network_path, paste0("cortab_", nm, ".csv"))
  write.csv(
    mat,
    file = out_file,
    row.names = TRUE
  )
}



#-提取全部图片的存储对象
plot = tab.r[[1]]
# 2--保存网络图可视化结果#-------
p0 = plot[[1]]
# 保存网络分析主图
ggsave(file.path(amplicon_network_path, "network_main.pdf"), plot = p0, width = 18, height = 6)


# 3--保存分组均值数据#------
tax = psg %>% vegan_tax() %>% as.data.frame()
head(tax)
tax$ID = row.names(tax)
dat = psg %>% ps.group.mean()
head(dat)
dat2 = dat %>% left_join(tax,by = "ID")
write_excel_csv(dat2,file.path(amplicon_network_path, "ps.group.mean.tax.csv"))

#4--网络属性#---------
cor = cortab
id = names(cor)
id
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

dat2
dat2 = as.data.frame(dat2)
dat2$ID = row.names(dat2)
write_excel_csv(dat2,paste0(amplicon_network_path,"/net.propertits.csv"))


#5.0 --共享边#------
mats <- cortab
# ab <- align_corr_pair(mats$G, mats$B,
#                       mode = "union", fill = NA, name_std = F)
ab =   align_corr_pair(mats[[1]], mats[[2]], name_std = F)
# 例：|r|≥0.5 视为一条边；共享边不要求同号
res_net <- compare_corr_two(ab[[1]], ab[[2]], edge_thr = 0.6, edge_rule = "abs", same_sign_only = FALSE)

res_net$summary[c("A_edge","B_edge","shared_edges","A_only_edges","B_only_edges")]
head(res_net$shared_edges)


#5.1 统计共享边参与的节点及其占据的丰度比例-----
asv_vec = c(res_net$shared_edges$var1,res_net$shared_edges$var2) %>% unique()

tax_bac <- phyloseq::tax_table(psg) %>% as.data.frame()

# 2. 统计共享节点 & 门/属构成
res_shared <- get_shared_nodes_from_align(
  align_res  = ab,
  tax_table  = tax_bac,
  id_col     = NULL,
  phylum_col = "Phylum",
  genus_col  = "Genus"
)

# 共享网络里的所有节点
res_shared$shared_nodes

# 共享网络中匹配到 taxonomy 的微生物节点
res_shared$node_tax

# 共享节点主要来自哪些门？
res_shared$phylum_summary

# 主要是哪些属？
res_shared$genus_summary

#--5.2 统计共享特有节点分类信息#-----

tax_bac <- phyloseq::tax_table(psg) %>% as.data.frame()
head(tax_bac,100)
tax_bac$geneid = row.names(tax_bac)


tax_bac2 <- tax_bac %>%
  mutate(
    Phylum = Phylum,
    Genus  = NA_character_
  )

comp3 <- compare_shared_A_B_tax(
  res_net,
  tax_table  = tax_bac,
  id_col     = "geneid",
  phylum_col = "Phylum",
  genus_col  = "Class"
)


# 已确认，没问题
comp3$shared$nodes %>% length()
comp3$A_only$nodes%>% length()
comp3$B_only$nodes%>% length()

c(res_net$WT_only_edges$var1,res_net$WT_only_edges$var2) %>% unique() %>% length()



# 三类边中，参与节点来自哪些门（各自占比）
comp3$shared$phylum_summary
comp3$A_only$phylum_summary
comp3$B_only$phylum_summary

# 同一个门在 shared / A_only / B_only 中的节点数 & 相对比例，对比表：
comp3$phylum_compare %>% head(20)

# 同一个属在三类网络中的情况：
comp3$genus_compare %>% head(20)


tax_bac <- phyloseq::tax_table(psg) %>% as.data.frame()
tax_bac$geneid = row.names(tax_bac)
head(tax_bac)

#--5.3共享特有边来自跨门还是同门#--------

comp_phylum <- compare_phylum_structure_three(
  res_net    = res_net,
  tax_table  = tax_bac,
  var1_col   = "var1",
  var2_col   = "var2",
  id_col     = "geneid",   # 如果 ASV ID 在 geneid 列
  phylum_col = "Class2"
)

comp_phylum <- compare_phylum_structure_three(
  res_net,
  tax_table  = tax_bac,
  id_col     = "geneid",
  phylum_col = "Phylum"
)

# 1）看三类边整体上“同门 vs 跨门”的比例
comp_phylum$within_between_compare

# 2）shared 网络里，同门互作主要来自哪些门：
comp_phylum$shared$within_phylum

#   跨门互作主要是哪些门–门：
comp_phylum$shared$between_pairs

# 3）A_only（健康特异互作）同门/跨门的结构：
comp_phylum$A_only$within_phylum
comp_phylum$A_only$between_pairs

# 4）B_only（发病特异互作）同门/跨门的结构：
comp_phylum$B_only$within_phylum
comp_phylum$B_only$between_pairs



# 6--参与网络的共享特有节点#-----
nodes_A <- get_network_nodes(ab[[1]])
nodes_B <- get_network_nodes(ab[[2]])
# 各自网络中参与的节点（去掉相关性为 0 和 NA 后）
length(nodes_A); head(nodes_A)
length(nodes_B); head(nodes_B)

# 共同参与的节点（两个网络里都有边的）
shared_nodes <- intersect(nodes_A, nodes_B)
shared_nodes
# 仅在 A 网络参与的节点
A_only_nodes <- setdiff(nodes_A, nodes_B)
A_only_nodes
# 仅在 B 网络参与的节点
B_only_nodes <- setdiff(nodes_B, nodes_A)
B_only_nodes
# 如果你想看这三类里面各自多少个：
length(shared_nodes)
length(A_only_nodes)
length(B_only_nodes)



#7--韦恩图展示相同和不同的边#---------

library(VennDiagram); library(grid)

S <- res_net$summary
A  <- unname(S$A_edge)
B  <- unname(S$B_edge)
AB <- unname(S$shared_edges)

grid.newpage()
venn_grobs <-draw.pairwise.venn(
  area1 = A,
  area2 = B,
  cross.area = AB,
  category = c("Network A", "Network B"),
  fill = c("#4C97FF", "#FF8A5B"),
  alpha = 0.6,
  lwd = 1.2,
  cex = 1.2,
  cat.cex = 1.2
)
p_venn <- ggplotify::as.ggplot(grid::grobTree(venn_grobs))

## 4) 保存到与主图相同目录（文件名自定）
ggsave(
  filename = file.path(amplicon_network_path, "network_venn.pdf"),
  plot     = p_venn,
  width    = 6, height = 6
)



#8.1--共享网络展示#-------
tax_df <- as.data.frame(phyloseq::tax_table(psg))
tax_df$ASV <- rownames(tax_df)

# 保留常见层级
keep_ranks <- intersect(c("Kingdom","Phylum","Class","Order","Family","Genus"), colnames(tax_df))
tax_df <- tax_df[, c("ASV", keep_ranks), drop = FALSE]

# 共享网络可视化
edges <- res_net$shared_edges

# 左连接注释
edges_ann <- edges %>%
  left_join(tax_df, by = c("var1" = "ASV")) %>%
  rename_with(~ paste0(., "_1"), keep_ranks) %>%
  left_join(tax_df, by = c("var2" = "ASV")) %>%
  rename_with(~ paste0(., "_2"), keep_ranks)
head(edges_ann)
# 必要包
library(ggraph)
library(tidygraph)

ed_simple <- edges_ann %>%
  transmute(a = var1, b = var2, r = r_WT) %>%
  group_by(a,b) %>% summarise(r = mean(r, na.rm = TRUE), .groups = "drop")

nodes <- sort(unique(c(ed_simple$a, ed_simple$b)))
corA <- matrix(0, nrow = length(nodes), ncol = length(nodes),
               dimnames = list(nodes, nodes))
# 填充上三角
idx <- cbind(match(ed_simple$a, nodes), match(ed_simple$b, nodes))
corA[idx] <- ed_simple$r
# 对称化
corA <- corA + t(corA)
diag(corA) <- 0
dim(corA)

result2 = model_maptree2(cor = corA,
                         method = "cluster_fast_greedy"
)

node = result2[[1]]
head(node)
# # ---node节点注释

otu_table = psg %>% vegan_otu() %>%  t() %>% as.data.frame()
tax = psg %>% vegan_tax() %>% as.data.frame()
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax)
head(nodes)
dim(nodes)
nodes$elements %>% unique()
#-----计算边
edge = edgeBuild(cor = corA,node = node)

### 出图
pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                                data = edge, size = 1,alpha = 0.1) +
  geom_point(aes(X1, X2,fill = Phylum,size = mean),pch = 21, data = nodes) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

pnet

ggsave(file.path(amplicon_network_path, "network_共享网络图.pdf"), plot = pnet, width = 12, height = 7)


#8.2--共享网络门-class-属水平展示#----------

genus_freq <- c(edges_ann$Phylum_1, edges_ann$Phylum_1)
genus_freq <- table(genus_freq)
genus_freq <- sort(genus_freq, decreasing = TRUE)
head(genus_freq, 20)
library(ggplot2)

top_genus <- as.data.frame(genus_freq)
colnames(top_genus) <- c("Genus","Edges")

p = ggplot(top_genus, aes(x = reorder(Genus, Edges), y = Edges)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top genera involved in shared edges",
       x = "Genus", y = "Number of shared edges")

p
ggsave(file.path(amplicon_network_path, "network_共享网络图-组成门水平.pdf"), plot = p, width = 6, height = 6)


library(ggplot2)
genus_freq <- c(edges_ann$Genus_1, edges_ann$Genus_2)
genus_freq <- table(genus_freq)
genus_freq <- sort(genus_freq, decreasing = TRUE)
head(genus_freq, 20)

top_genus <- as.data.frame(genus_freq)
colnames(top_genus) <- c("Genus","Edges")

p = ggplot(top_genus, aes(x = reorder(Genus, Edges), y = Edges)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top genera involved in shared edges",
       x = "Genus", y = "Number of shared edges")

ggsave(file.path(amplicon_network_path, "network_共享网络图-组成属水平.pdf"), plot = p, width = 6, height = 6)



library(ggplot2)
# genus_freq <- c(edges_ann$Order_1, edges_ann$Order_2)
genus_freq <- c(edges_ann$Class_1, edges_ann$Class_1)


genus_freq <- table(genus_freq)
genus_freq <- sort(genus_freq, decreasing = TRUE)
head(genus_freq, 20)

top_genus <- as.data.frame(genus_freq)
colnames(top_genus) <- c("Genus","Edges")

p = ggplot(top_genus, aes(x = reorder(Genus, Edges), y = Edges)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top genera involved in shared edges",
       x = "Genus", y = "Number of shared edges")
p
ggsave(file.path(amplicon_network_path, "network_共享网络图-组成class水平.pdf"), plot = p, width = 6, height = 6)




#8.3--共享网络微生物堆叠柱状图#----

tab =  res_net$shared_edges
id = c(tab$var1,tab$var2) %>% unique()

#8.4--barMainplot.micro: 堆积柱状图展示组成----
library(ggalluvial)
pst = ps.16s %>%
  scale_micro() %>%
  subset_taxa.wt("OTU",id)

result = barMainplot.micro(ps = pst,
                           j = "Genus",
                           # axis_ord = axis_order,
                           label = FALSE,
                           sd = FALSE,
                           tran = F,
                           Top =10)

p4_1 = result[[1]] +
  scale_fill_manual(values = colset2) +
  scale_x_discrete(limits = axis_order) +theme_nature()+
  theme(axis.title.y = element_text(angle = 90))

p4_1

ggsave(file.path(amplicon_network_path, "network_共享网络图-属水平堆叠柱状图.pdf"), plot = p4_1 , width = 6, height = 8)
result = barMainplot.micro(ps = pst,
                           j = "Phylum",
                           # axis_ord = axis_order,
                           label = FALSE,
                           sd = FALSE,
                           tran = F,
                           Top =10)

p4_1 = result[[1]] +
  scale_fill_manual(values = colset2) +
  scale_x_discrete(limits = axis_order) +theme_nature()+
  theme(axis.title.y = element_text(angle = 90))

ggsave(file.path(amplicon_network_path, "network_共享网络图-门水平堆叠柱状图.pdf"), plot = p4_1 , width = 6, height = 8)


#9.1一组中特有的网络关系#------

edges <- res_net$WT_only_edges

# 左连接注释
edges_ann <- edges %>%
  left_join(tax_df, by = c("var1" = "ASV")) %>%
  rename_with(~ paste0(., "_1"), keep_ranks) %>%
  left_join(tax_df, by = c("var2" = "ASV")) %>%
  rename_with(~ paste0(., "_2"), keep_ranks)

head(edges_ann)
dim(edges_ann)

ed_simple <- edges_ann %>%
  transmute(a = var1, b = var2, r = r_WT) %>%
  group_by(a,b) %>% summarise(r = mean(r, na.rm = TRUE), .groups = "drop")

nodes <- sort(unique(c(ed_simple$a, ed_simple$b)))
corA <- matrix(0, nrow = length(nodes), ncol = length(nodes),
               dimnames = list(nodes, nodes))
# 填充上三角
idx <- cbind(match(ed_simple$a, nodes), match(ed_simple$b, nodes))
corA[idx] <- ed_simple$r
# 对称化
corA <- corA + t(corA)
diag(corA) <- 0

dim(corA)

result2 = model_maptree2(cor = corA,
                         method = "cluster_fast_greedy"
)

# result2 = PolygonRrClusterG (cor = cor,nodeGroup =group2 )
node = result2[[1]]
head(node)
# # ---node节点注释

otu_table = psg %>% vegan_otu() %>%  t() %>% as.data.frame()
tax = psg %>% vegan_tax() %>% as.data.frame()
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax)
head(nodes)
#
# nodes2 = nodes %>% inner_join(netClu,by = c("elements" = "ID"))
# nodes2$group = paste("Model_",nodes2$group,sep = "")

#-----计算边
edge = edgeBuild(cor = corA,node = node)

### 出图
pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                                data = edge, size = 1,alpha = 0.1) +
  geom_point(aes(X1, X2,fill = Phylum,size = mean),pch = 21, data = nodes) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

pnet

ggsave(file.path(amplicon_network_path, "network_组1.pdf"), plot = pnet , width = 15, height = 7)


# 共享网络主要有哪些微生物构成

genus_freq <- c(edges_ann$Phylum_1, edges_ann$Phylum_1)
genus_freq <- table(genus_freq)
genus_freq <- sort(genus_freq, decreasing = TRUE)
head(genus_freq, 20)


top_genus <- as.data.frame(genus_freq)
colnames(top_genus) <- c("Genus","Edges")

p = ggplot(top_genus, aes(x = reorder(Genus, Edges), y = Edges)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top genera involved in shared edges",
       x = "Genus", y = "Number of shared edges")
p
ggsave(file.path(amplicon_network_path, "network_健康特有网络图-门水平组成.pdf"), plot = p , width = 6, height = 6)


library(ggplot2)
genus_freq <- c(edges_ann$Genus_1, edges_ann$Genus_2)
genus_freq <- table(genus_freq)
genus_freq <- sort(genus_freq, decreasing = TRUE)
head(genus_freq, 20)

top_genus <- as.data.frame(genus_freq)[1:20,]
colnames(top_genus) <- c("Genus","Edges")

p  = ggplot(top_genus, aes(x = reorder(Genus, Edges), y = Edges)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top genera involved in shared edges",
       x = "Genus", y = "Number of shared edges")

p

ggsave(file.path(amplicon_network_path, "network_健康特有网络图-属水平组成.pdf"), plot = p, width = 6, height = 6)

genus_freq <- c(edges_ann$Class_1, edges_ann$Class_2)
genus_freq <- table(genus_freq)
genus_freq <- sort(genus_freq, decreasing = TRUE)
head(genus_freq, 20)

top_genus <- as.data.frame(genus_freq)[1:15,]
colnames(top_genus) <- c("Genus","Edges")

p  = ggplot(top_genus, aes(x = reorder(Genus, Edges), y = Edges)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top genera involved in shared edges",
       x = "Genus", y = "Number of shared edges")

p


ggsave(file.path(amplicon_network_path, "network_健康特有网络图-class水平组成.pdf"), plot = p , width = 6, height = 6)




# 2 共享网络的微生物组成是什么情况

tab =  res_net$WT_only_edges
id = c(tab$var1,tab$var2) %>% unique()

#16 barMainplot.micro: 堆积柱状图展示组成

library(ggalluvial)
pst = ps.16s %>%
  scale_micro() %>%
  subset_taxa.wt("OTU",id)

result = barMainplot.micro(ps = pst,
                           j = "Phylum",
                           # axis_ord = axis_order,
                           label = FALSE,
                           sd = FALSE,
                           tran = F,
                           Top =10)
p4_1 <- result[[1]]+
  scale_fill_manual(values = colset2) +
  scale_x_discrete(limits = axis_order) +theme_nature()+
  theme(axis.title.y = element_text(angle = 90))

p4_1
ggsave(file.path(amplicon_network_path, "network_健康特有网络图-门水平堆叠柱状图.pdf"), plot = p4_1 , width = 6, height = 8)



result = barMainplot.micro(ps = pst,
                           j = "Genus",
                           # axis_ord = axis_order,
                           label = FALSE,
                           sd = FALSE,
                           tran = F,
                           Top =10)
p4_1 <- result[[1]]+
  scale_fill_manual(values = colset2) +
  scale_x_discrete(limits = axis_order) +theme_nature()+
  theme(axis.title.y = element_text(angle = 90))


ggsave(file.path(amplicon_network_path, "network_健康特有网络图-属水平堆叠柱状图.pdf"), plot = p4_1 , width = 6, height = 8)




result = barMainplot.micro(ps = pst,
                           j = "Class",
                           # axis_ord = axis_order,
                           label = FALSE,
                           sd = FALSE,
                           tran = F,
                           Top =10)
p4_1 <- result[[1]]+
  scale_fill_manual(values = colset2) +
  scale_x_discrete(limits = axis_order) +theme_nature()+
  theme(axis.title.y = element_text(angle = 90))


ggsave(file.path(amplicon_network_path, "network_健康特有网络图-class水平堆叠柱状图.pdf"), plot = p4_1 , width = 6, height = 8)

result = barMainplot.micro(ps = pst,
                           j = "Family",
                           # axis_ord = axis_order,
                           label = FALSE,
                           sd = FALSE,
                           tran = F,
                           Top =11)
p4_1 <- result[[1]]+
  scale_fill_manual(values = colset2) +
  scale_x_discrete(limits = axis_order) +theme_nature()+
  theme(axis.title.y = element_text(angle = 90))


#10.1 另一组网络中特有的网络关系#----

edges <- res_net$OE_only_edges

# 左连接注释
edges_ann <- edges %>%
  left_join(tax_df, by = c("var1" = "ASV")) %>%
  rename_with(~ paste0(., "_1"), keep_ranks) %>%
  left_join(tax_df, by = c("var2" = "ASV")) %>%
  rename_with(~ paste0(., "_2"), keep_ranks)

head(edges_ann)
dim(edges_ann)

ed_simple <- edges_ann %>%
  transmute(a = var1, b = var2, r = r_OE) %>%
  group_by(a,b) %>% summarise(r = mean(r, na.rm = TRUE), .groups = "drop")

nodes <- sort(unique(c(ed_simple$a, ed_simple$b)))
corA <- matrix(0, nrow = length(nodes), ncol = length(nodes),
               dimnames = list(nodes, nodes))
# 填充上三角
idx <- cbind(match(ed_simple$a, nodes), match(ed_simple$b, nodes))
corA[idx] <- ed_simple$r
# 对称化
corA <- corA + t(corA)
diag(corA) <- 0

dim(corA)

result2 = model_maptree2(cor = corA,
                         method = "cluster_fast_greedy"
)

node = result2[[1]]
head(node)


otu_table = psg %>% vegan_otu() %>%  t() %>% as.data.frame()
tax = psg %>% vegan_tax() %>% as.data.frame()
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax)
head(nodes)

#-----计算边
edge = edgeBuild(cor = corA,node = node)

### 出图
pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                                data = edge, size = 1,alpha = 0.1) +
  geom_point(aes(X1, X2,fill = Phylum,size = mean),pch = 21, data = nodes) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

pnet


ggsave(file.path(amplicon_network_path, "network_发病特有网络图.pdf"), plot = pnet , width = 15, height = 7)


# 共享网络主要有哪些微生物构成

genus_freq <- c(edges_ann$Phylum_1, edges_ann$Phylum_1)
genus_freq <- table(genus_freq)
genus_freq <- sort(genus_freq, decreasing = TRUE)
head(genus_freq, 20)
library(ggplot2)

top_genus <- as.data.frame(genus_freq)[1:15,]
colnames(top_genus) <- c("Genus","Edges")

p = ggplot(top_genus, aes(x = reorder(Genus, Edges), y = Edges)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top genera involved in shared edges",
       x = "Genus", y = "Number of shared edges")

p
ggsave(file.path(amplicon_network_path, "network_发病特有网络图-门水平组成.pdf"), plot = p , width = 6, height = 6)


genus_freq <- c(edges_ann$Genus_1, edges_ann$Genus_2)
genus_freq <- table(genus_freq)
genus_freq <- sort(genus_freq, decreasing = TRUE)
head(genus_freq, 20)

top_genus <- as.data.frame(genus_freq)[1:15,]
colnames(top_genus) <- c("Genus","Edges")

p = ggplot(top_genus, aes(x = reorder(Genus, Edges), y = Edges)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top genera involved in shared edges",
       x = "Genus", y = "Number of shared edges")
p

ggsave(file.path(amplicon_network_path, "network_发病特有网络图-属水平组成.pdf"), plot = p , width = 6, height = 6)


genus_freq <- c(edges_ann$Class_1, edges_ann$Class_2)
genus_freq <- table(genus_freq)
genus_freq <- sort(genus_freq, decreasing = TRUE)
head(genus_freq, 20)

top_genus <- as.data.frame(genus_freq)[1:15,]
colnames(top_genus) <- c("Genus","Edges")

p = ggplot(top_genus, aes(x = reorder(Genus, Edges), y = Edges)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top genera involved in shared edges",
       x = "Genus", y = "Number of shared edges")
p

ggsave(file.path(amplicon_network_path, "network_发病特有网络图-class水平组成.pdf"), plot = p , width = 6, height = 6)


# 2 共享网络的微生物组成是什么情况

tab =  res_net$OE_only_edges
id = c(tab$var1,tab$var2) %>% unique()


pst = ps.16s %>%
  scale_micro() %>%
  subset_taxa.wt("OTU",id)

result = barMainplot.micro(ps = pst,
                           j = "Phylum",
                           # axis_ord = axis_order,
                           label = FALSE,
                           sd = FALSE,
                           tran = F,
                           Top =10)
p4_1 <- result[[1]]+
  scale_fill_manual(values = colset2) +
  scale_x_discrete(limits = axis_order) +theme_nature()+
  theme(axis.title.y = element_text(angle = 90))


ggsave(file.path(amplicon_network_path, "network_发病特有网络图-门水平堆叠柱状图.pdf"), plot = p4_1 , width = 6, height = 8)

result = barMainplot.micro(ps = pst,
                           j = "Genus",
                           # axis_ord = axis_order,
                           label = FALSE,
                           sd = FALSE,
                           tran = F,
                           Top =10)
p4_1 <- result[[1]]+
  scale_fill_manual(values = colset2) +
  scale_x_discrete(limits = axis_order) +theme_nature()+
  theme(axis.title.y = element_text(angle = 90))

p4_1

ggsave(file.path(amplicon_network_path, "network_发病特有网络图-属水平堆叠柱状图.pdf"), plot = p4_1 , width = 6, height = 8)






