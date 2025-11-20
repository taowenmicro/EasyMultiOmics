
find_monotonic_taxa <- function(ps, group = "Group", method = c("spearman","strict")) {
  method <- match.arg(method)

  # 1. 提取 OTU 表
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu <- t(otu)

  # 2. 提取分组信息
  meta <- as.data.frame(sample_data(ps))
  groups <- meta[[group]]

  # 3. 按分组求均值
  otu_sum <- rowsum(otu, groups)
  otu_mean <- otu_sum / as.vector(table(groups))
  otu_mean_df <- data.frame(Group = rownames(otu_mean), otu_mean, check.names = FALSE)

  # 分组编号 (1,2,3,4,5...) —— 必须是有序的
  group_levels <- seq_len(nrow(otu_mean_df))

  # 只保留数值矩阵
  otu_mat <- as.matrix(otu_mean_df[,-1])

  # 4. 初始化结果
  increasing <- decreasing <- character()

  # 5. 遍历每个 OTU
  for (i in seq_len(ncol(otu_mat))) {
    vals <- otu_mat[, i]
    if (all(is.na(vals)) || var(vals, na.rm = TRUE) == 0) next

    if (method == "spearman") {
      rho <- suppressWarnings(cor.test(group_levels, vals, method = "spearman", exact = FALSE))
      if (!is.na(rho$p.value) && rho$p.value < 0.05 && rho$estimate > 0) {
        increasing <- c(increasing, colnames(otu_mat)[i])
      }
      if (!is.na(rho$p.value) && rho$p.value < 0.05 && rho$estimate < 0) {
        decreasing <- c(decreasing, colnames(otu_mat)[i])
      }
    }

    if (method == "strict") {
      if (all(diff(vals) > 0)) increasing <- c(increasing, colnames(otu_mat)[i])
      if (all(diff(vals) < 0)) decreasing <- c(decreasing, colnames(otu_mat)[i])
    }
  }

  # 6. 结果表格：提取增加/减少 OTU 的丰度均值
  inc_df <- otu_mean_df %>% dplyr::select(Group, all_of(increasing))
  dec_df <- otu_mean_df %>% dplyr::select(Group, all_of(decreasing))
  head(inc_df)
  row.names(inc_df) = inc_df$Group
  inc_df$Group = NULL
  inc_df1 = t(inc_df) %>% as.data.frame()

  row.names(dec_df) = dec_df$Group
  dec_df$Group = NULL
  dec_df1 = t(dec_df) %>% as.data.frame()

  return(list(
    increasing_taxa = increasing,
    decreasing_taxa = decreasing,
    increasing_abundance = inc_df1,
    decreasing_abundance = dec_df1
  ))
}
