#' @title  Retrieve metabolite information from the HMDB database
#' @description
#' The ann.HMDB function retrieves relevant information about a metabolite based
#'  on its ID from the HMDB (Human Metabolome Database) stored in the db.metabolites database.
#' @param id The metabolite ID for which information is to be retrieved.
#' @return A data frame containing the metabolite information including HMDB ID, KEGG ID,
#' Super_class, Class, Sub_class and related details from the database.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njau.edu.cn}
#' @examples
#' library(dplyr)
#' tax= ps.ms %>% vegan_tax() %>%as.data.frame()
#' id = tax$Metabolite
#' tax1 = ann.HMDB (id = id)
#' head(tax1)
#' @export
ann.HMDB <- function(id, verbose = FALSE) {

  db <- get("db.metabolites", envir = .GlobalEnv)

  # 预处理数据库：创建小写版本和清洗版本
  if(!"Name_Lower" %in% colnames(db)) {
    db$Name_Lower <- str_to_lower(db$Name)
  }

  if(!"Name_Clean" %in% colnames(db)) {
    db$Name_Clean <- clean_metabolite_name(db$Name)
  }

  # 初始化结果
  res_df <- data.frame(id.org = id, stringsAsFactors = FALSE)
  match_indices <- rep(NA, length(id))
  match_method <- rep(NA_character_, length(id))

  # 精确匹配（原始名称）
  id_lower <- str_to_lower(str_trim(id))
  idx_exact <- match(id_lower, db$Name_Lower)

  matched <- !is.na(idx_exact)
  match_indices[matched] <- idx_exact[matched]
  match_method[matched] <- "exact"

  if(verbose) cat("精确匹配:", sum(matched), "/", length(id), "\n")

  # 去除末尾数字
  na_locs <- which(is.na(match_indices))
  if(length(na_locs) > 0) {
    id_no_num <- id[na_locs] %>%
      str_replace_all("[ ][0-9]+$", "") %>%
      str_to_lower() %>%
      str_trim()

    idx_no_num <- match(id_no_num, db$Name_Lower)
    matched <- !is.na(idx_no_num)
    match_indices[na_locs[matched]] <- idx_no_num[matched]
    match_method[na_locs[matched]] <- "remove_trailing_num"

    if(verbose) cat("去除末尾数字后匹配:", sum(matched), "个\n")
  }

  # 去除立体化学前缀
  na_locs <- which(is.na(match_indices))
  if(length(na_locs) > 0) {
    id_no_stereo <- id[na_locs] %>%
      str_replace_all(regex("^(l|d|alpha|beta|cis|trans|n|o|s|p|m)-", ignore_case = TRUE), "") %>%
      str_replace_all(regex("^(dl|rac)-", ignore_case = TRUE), "") %>%
      str_to_lower() %>%
      str_trim()

    idx_no_stereo <- match(id_no_stereo, db$Name_Lower)
    matched <- !is.na(idx_no_stereo)
    match_indices[na_locs[matched]] <- idx_no_stereo[matched]
    match_method[na_locs[matched]] <- "remove_stereo"

    if(verbose) cat("去除立体化学前缀后匹配:", sum(matched), "个\n")
  }

  # ate <-> ic acid 转换
  na_locs <- which(is.na(match_indices))
  if(length(na_locs) > 0) {
    id_acid_conv <- id[na_locs] %>%
      str_replace_all(regex("^(l|d|alpha|beta|cis|trans|n|o|s|p|m|dl|rac)-", ignore_case = TRUE), "") %>%
      str_to_lower() %>%
      str_replace_all("ate$", "ic acid") %>%
      str_trim()

    idx_acid <- match(id_acid_conv, db$Name_Lower)
    matched <- !is.na(idx_acid)
    match_indices[na_locs[matched]] <- idx_acid[matched]
    match_method[na_locs[matched]] <- "ate_to_acid"

    if(verbose) cat("ate转ic acid后匹配:", sum(matched), "个\n")
  }

  # ic acid <-> ate 反向转换
  na_locs <- which(is.na(match_indices))
  if(length(na_locs) > 0) {
    id_ate_conv <- id[na_locs] %>%
      str_replace_all(regex("^(l|d|alpha|beta|cis|trans|n|o|s|p|m|dl|rac)-", ignore_case = TRUE), "") %>%
      str_to_lower() %>%
      str_replace_all("ic acid$", "ate") %>%
      str_trim()

    idx_ate <- match(id_ate_conv, db$Name_Lower)
    matched <- !is.na(idx_ate)
    match_indices[na_locs[matched]] <- idx_ate[matched]
    match_method[na_locs[matched]] <- "acid_to_ate"

    if(verbose) cat("ic acid转ate后匹配:", sum(matched), "个\n")
  }

  # 深度清洗匹配
  na_locs <- which(is.na(match_indices))
  if(length(na_locs) > 0) {
    id_deep_clean <- clean_metabolite_name(id[na_locs])

    idx_clean <- match(id_deep_clean, db$Name_Clean)
    matched <- !is.na(idx_clean)
    match_indices[na_locs[matched]] <- idx_clean[matched]
    match_method[na_locs[matched]] <- "deep_clean"

    if(verbose) cat("深度清洗后匹配:", sum(matched), "个\n")
  }

  # 去除括号内容
  na_locs <- which(is.na(match_indices))
  if(length(na_locs) > 0) {
    id_no_paren <- id[na_locs] %>%
      str_replace_all("\\([^)]*\\)", "") %>%
      str_replace_all("\\[[^]]*\\]", "") %>%
      str_to_lower() %>%
      str_trim()

    idx_no_paren <- match(id_no_paren, db$Name_Lower)
    matched <- !is.na(idx_no_paren)
    match_indices[na_locs[matched]] <- idx_no_paren[matched]
    match_method[na_locs[matched]] <- "remove_parentheses"

    if(verbose) cat("去除括号后匹配:", sum(matched), "个\n")
  }

  # 模糊匹配 (agrep)
  na_locs <- which(is.na(match_indices))
  if(length(na_locs) > 0 && length(na_locs) <= 100) {  # 限制数量避免太慢
    for(i in na_locs) {
      query <- str_to_lower(str_trim(id[i]))
      if(nchar(query) < 4) next  # 太短的名称不做模糊匹配

      # 使用 agrep 进行模糊匹配，最大编辑距离为2
      fuzzy_match <- agrep(query, db$Name_Lower, max.distance = 2, value = FALSE)

      if(length(fuzzy_match) == 1) {  # 只有唯一匹配时才接受
        match_indices[i] <- fuzzy_match[1]
        match_method[i] <- "fuzzy"
      }
    }

    if(verbose) cat("模糊匹配:", sum(match_method == "fuzzy", na.rm = TRUE), "个\n")
  }

  # 提取匹配结果
  # 处理匹配结果，保持与输入相同的行数
  db_matched <- db[match_indices, , drop = FALSE]

  # 对于未匹配的行，填充NA
  if(any(is.na(match_indices))) {
    # 确保维度正确
    na_rows <- which(is.na(match_indices))
    for(col in names(db)) {
      if(!col %in% names(db_matched)) {
        db_matched[[col]] <- rep(NA, nrow(db_matched))
      }
    }
  }

  db_matched$match_method <- match_method

  # 合并结果
  final_res <- cbind(res_df, db_matched)

  # 统计
  if(verbose) {
    cat("\n========== 匹配统计 ==========\n")
    cat("总数:", length(id), "\n")
    cat("匹配成功:", sum(!is.na(match_indices)), "\n")
    cat("匹配率:", round(sum(!is.na(match_indices))/length(id)*100, 2), "%\n")
    cat("\n各策略匹配数:\n")
    print(table(match_method, useNA = "ifany"))
  }

  return(final_res)
}

# 辅助函数：深度清洗代谢物名称
clean_metabolite_name <- function(name) {
  name %>%
    as.character() %>%
    # 去除立体化学信息
    str_replace_all(regex("^(l|d|dl|rac|alpha|beta|cis|trans|n|o|s|p|m)-", ignore_case = TRUE), "") %>%
    str_replace_all(regex("\\((l|d|dl|rac|alpha|beta|cis|trans|r|s|e|z)\\)", ignore_case = TRUE), "") %>%
    # 去除数字编号
    str_replace_all("[ ][0-9]+$", "") %>%
    # 标准化分隔符
    str_replace_all("[_-]+", " ") %>%
    # 去除多余空格
    str_replace_all("\\s+", " ") %>%
    # 转小写并去除首尾空格
    str_to_lower() %>%
    str_trim()
}
