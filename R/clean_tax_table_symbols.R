
# library(phyloseq)
#
# # 假设你有一个 ps.micro
# ps.cs <- clean_tax_table_symbols(ps.cs)
#
# # 看看修改后的 tax 表
# head(as.data.frame(phyloseq::tax_table(ps.micro.clean)))

clean_tax_table_symbols <- function(ps,
                                    replacement = ".",
                                    allowed_chars = "A-Za-z0-9_.") {
  if (!inherits(ps, "phyloseq")) {
    stop("ps must be a phyloseq object.")
  }

  # 取出 tax_table
  tax <- phyloseq::tax_table(ps)

  # 转成 matrix，方便逐元素处理
  tax_mat <- as.matrix(tax)

  # 清洗函数：处理一个字符向量
  clean_fun <- function(x) {
    x_chr <- as.character(x)

    # 保留 NA
    na_idx <- is.na(x_chr)

    # 1) 先把 "-" 替换为 "."
    x_chr[!na_idx] <- gsub("-", replacement, x_chr[!na_idx], fixed = TRUE)

    # 2) 再把除 allowed_chars 外的所有字符替换为 "."
    pattern <- sprintf("[^%s]", allowed_chars)
    x_chr[!na_idx] <- gsub(pattern, replacement, x_chr[!na_idx])

    x_chr
  }

  # 对 tax 表每一列进行清洗
  tax_clean <- apply(tax_mat, 2, clean_fun)

  # 保持原来的 dimnames
  rownames(tax_clean) <- rownames(tax_mat)
  colnames(tax_clean) <- colnames(tax_mat)

  # 写回 phyloseq 对象
  phyloseq::tax_table(ps) <- phyloseq::tax_table(tax_clean)

  return(ps)
}
