# 辅助：根据 kingdom/domain 字符串判断属于哪个界 -------------------
.classify_kingdom_to_type <- function(x) {
  x <- gsub("^[a-z]__", "", x)   # 去掉 k__ 这类前缀
  x <- trimws(x)
  
  res <- rep(NA_character_, length(x))
  
  # 细菌 / 古菌
  res[grepl("bacteria|archaea", x, ignore.case = TRUE)] <- "bacteria"
  
  # 真菌
  res[grepl("fungi|ascomycota|basidiomycota", x, ignore.case = TRUE)] <- "fungi"
  
  # 病毒
  res[grepl("virus|viridae|phage", x, ignore.case = TRUE)] <- "virus"
  
  # 原生动物 / 原生生物
  res[grepl("protist|protozoa|ciliophora|apicomplexa|amoebozoa",
            x, ignore.case = TRUE)] <- "protist"
  
  res
}

# 计算不同界的丰度比例（按 reads/序列数） -----------------------------
# 返回一个 named integer 向量，如：
#   bacteria fungi
#        60    40
get_amplicon_type_ratio <- function(ps) {
  # 1. 取 OTU 表
  otu <- phyloseq::otu_table(ps)
  otu_mat <- as(otu, "matrix")
  
  # 行是否为 taxon
  if (!phyloseq::taxa_are_rows(otu)) {
    otu_mat <- t(otu_mat)
  }
  
  # 2. 取 taxonomy
  tt <- phyloseq::tax_table(ps)
  if (is.null(tt)) {
    return(structure(100L, names = "amplicon"))
  }
  tt <- as.matrix(tt)
  
  col_lower <- tolower(colnames(tt))
  kingdom_cols <- which(col_lower %in% c("kingdom", "k", "domain", "superkingdom"))
  
  # 没有明确列名时，所有列都扫一遍
  if (length(kingdom_cols) == 0) {
    kingdom_cols <- seq_len(ncol(tt))
  }
  
  kingdom_raw <- apply(tt[, kingdom_cols, drop = FALSE], 1, function(v) {
    v <- v[!is.na(v)]
    if (length(v) == 0) return(NA_character_)
    v[1]
  })
  
  # 3. 对每个 taxon 判断所属界
  taxon_type <- .classify_kingdom_to_type(kingdom_raw)
  
  # 4. 计算每个界的总丰度
  type_order <- c("bacteria", "fungi", "protist", "virus")
  abund_vec  <- numeric(length(type_order))
  names(abund_vec) <- type_order
  
  for (tp in type_order) {
    sel <- which(taxon_type == tp)
    if (length(sel) > 0) {
      abund_vec[tp] <- sum(otu_mat[sel, , drop = FALSE])
    }
  }
  
  # 没有识别出任何界，就返回通用标签
  total_abund <- sum(abund_vec)
  if (total_abund <= 0) {
    return(structure(100L, names = "amplicon"))
  }
  
  # 5. 计算百分比（四舍五入），不带 % 符号
  pct <- round(abund_vec / total_abund * 100)
  
  # 去掉 0% 的界
  pct <- pct[pct > 0]
  
  pct
}

# 根据比例构建“高档”标签 ---------------------------------------------
# 单一界：  bacteria(100) -> "bacteria"
# 多界：    Bac60-Fun40-Pro20 -> "multi_Bac60-Fun40-Pro20"
build_amplicon_type_label_with_ratio <- function(type_pct) {
  # 单一类型直接返回英文名
  if (length(type_pct) == 1) {
    nm <- names(type_pct)
    if (identical(nm, "amplicon")) {
      return("amplicon_unclassified")
    }
    return(nm)
  }
  
  short_map <- c(
    bacteria = "Bac",
    fungi    = "Fun",
    protist  = "Pro",
    virus    = "Vir"
  )
  
  # 固定顺序，确保命名稳定
  type_order <- c("bacteria", "fungi", "protist", "virus")
  # 只保留我们识别的类型
  keep <- intersect(type_order, names(type_pct))
  type_pct <- type_pct[keep]
  
  if (length(type_pct) == 0) {
    return("amplicon_unclassified")
  }
  
  parts <- paste0(short_map[names(type_pct)], type_pct)
  
  paste0("multi_", paste(parts, collapse = "-"))
}

# 创建结果目录主函数（新版：时间可选） -------------------------------
create_amplicon_result_dir <- function(ps,
                                       base_dir = "./result",
                                       prefix   = "amplicon",
                                       include_time = FALSE) {
  # 1. 计算不同界在样本中的丰度比例
  type_pct <- get_amplicon_type_ratio(ps)
  
  # 2. 构建优雅的类型标签（单界 or 多界+比例）
  type_label <- build_amplicon_type_label_with_ratio(type_pct)
  # 例如：
  #   "bacteria"
  #   "multi_Bac60-Fun40"
  #   "multi_Bac50-Fun30-Pro20"
  
  # 3. 日期必带：年月日，例如 20251129
  date_part <- format(Sys.Date(), "%Y%m%d")
  
  # 4. 时间部分可选：时分秒，例如 103522
  if (isTRUE(include_time)) {
    time_part <- format(Sys.time(), "%H%M%S")
    datetime_label <- paste(date_part, time_part, sep = "_")
  } else {
    datetime_label <- date_part
  }
  
  # 5. 组合最终目录名
  #    不含时间：amplicon_multi_Bac60-Fun40_20251129
  #    含时间：  amplicon_multi_Bac60-Fun40_20251129_103522
  folder_name <- paste(prefix, type_label, datetime_label, sep = "_")
  out_path    <- file.path(base_dir, folder_name)
  
  # 6. 创建目录：优先用你已有的 dir_create()
  if (exists("dir_create")) {
    dir_create(out_path)
  } else if (requireNamespace("fs", quietly = TRUE)) {
    fs::dir_create(out_path, recurse = TRUE)
  } else {
    dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  }
  
  out_path
}
