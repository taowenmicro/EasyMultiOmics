#' @export
download_kegg_database <- function(verbose = TRUE) {

  if(verbose) cat("正在从KEGG REST API下载数据...\n")

  tryCatch({
    # 方法1: 使用KEGG REST API获取所有compound列表
    # 获取所有compound ID
    compound_list_url <- "https://rest.kegg.jp/list/compound"
    compound_data <- read.table(compound_list_url, sep = "\t",
                                quote = "", fill = TRUE,
                                stringsAsFactors = FALSE,
                                comment.char = "")

    if(verbose) cat("  获取到", nrow(compound_data), "个compounds\n")

    # 整理数据
    keggID <- compound_data[, 1]
    names_raw <- compound_data[, 2]

    # 解析名称
    metabolites_list <- strsplit(names_raw, ";\\s*")

    # 创建数据框
    max_names <- max(sapply(metabolites_list, length))

    # 初始化结果数据框
    result <- data.frame(
      keggID = keggID,
      allMetabolites = names_raw,
      stringsAsFactors = FALSE
    )

    # 添加分列的代谢物名称
    for(i in 1:max_names) {
      col_name <- paste0("metabolites", i)
      result[[col_name]] <- sapply(metabolites_list, function(x) {
        if(length(x) >= i) return(x[i]) else return("")
      })
    }

    # 添加数据库时间戳属性
    attr(result, "download_time") <- Sys.time()
    attr(result, "download_date") <- as.character(Sys.Date())
    attr(result, "source") <- "KEGG REST API"

    if(verbose) {
      cat("数据整理完成\n")
      cat("  总记录数:", nrow(result), "\n")
      cat("  总列数:", ncol(result), "\n")
      cat("  下载时间:", as.character(Sys.time()), "\n")
    }

    return(result)

  }, error = function(e) {
    if(verbose) cat("从KEGG REST API下载失败:", e$message, "\n")

    # 备用方法：如果全局环境有数据，使用它
    if(exists("db.ms.kegg", envir = .GlobalEnv)) {
      if(verbose) cat("使用全局环境中的 db.ms.kegg 作为备用\n")
      db <- get("db.ms.kegg", envir = .GlobalEnv)
      # 添加时间戳
      attr(db, "download_time") <- Sys.time()
      attr(db, "download_date") <- as.character(Sys.Date())
      attr(db, "source") <- "Global Environment (Backup)"
      return(db)
    }

    stop("无法下载KEGG数据库，且没有备用数据")
  })
}


# KEGG数据库下载和缓存管理函数
get_kegg_database <- function(cache_dir = "./kegg_cache",
                              max_age_days = 90,
                              force_update = FALSE,
                              verbose = TRUE) {

  # 创建缓存目录
  if(!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
    if(verbose) cat("创建缓存目录:", cache_dir, "\n")
  }

  # 检查现有缓存文件（任何日期的）
  existing_files <- list.files(cache_dir, pattern = "^kegg_database_\\d{8}\\.rds$", full.names = TRUE)

  need_update <- force_update
  db_cached <- NULL
  latest_file <- NULL

  if(!need_update && length(existing_files) > 0) {
    # 找到最新的文件
    file_dates <- gsub(".*_(\\d{8})\\.rds$", "\\1", basename(existing_files))
    latest_idx <- which.max(as.Date(file_dates, format = "%Y%m%d"))
    latest_file <- existing_files[latest_idx]
    latest_date <- file_dates[latest_idx]

    # 读取缓存
    db_cached <- readRDS(latest_file)

    # 计算文件年龄
    file_date <- as.Date(latest_date, format = "%Y%m%d")
    days_old <- as.numeric(Sys.Date() - file_date)

    if(verbose) {
      cat("========== KEGG数据库缓存状态 ==========\n")
      cat("发现本地缓存:", basename(latest_file), "\n")
      cat("缓存日期:", format(file_date, "%Y-%m-%d"), "\n")
    }

    if(days_old > max_age_days) {
      if(verbose) cat("缓存已过期（>", max_age_days, "天），需要更新\n\n")
      need_update <- TRUE
    } else {
      if(verbose) {
        cat("缓存仍有效（剩余", max_age_days - days_old, "天），使用本地数据\n")
        cat("========================================\n\n")
      }
      return(db_cached)
    }
  } else {
    if(verbose) {
      cat("========== KEGG数据库缓存状态 ==========\n")
      cat("未找到本地缓存，需要下载\n")
      cat("========================================\n\n")
    }
    need_update <- TRUE
  }

  # 如果需要更新，从在线下载
  if(need_update) {
    tryCatch({
      # 从KEGG REST API下载
      db <- download_kegg_database(verbose = verbose)

      # 生成带日期的文件名
      current_date <- format(Sys.Date(), "%Y%m%d")
      cache_file <- file.path(cache_dir, paste0("kegg_database_", current_date, ".rds"))

      # 保存到本地
      saveRDS(db, cache_file)

      if(verbose) {
        cat("\n========== 缓存保存成功 ==========\n")
        cat("数据库已保存到:", cache_file, "\n")
        cat("数据库大小:", format(object.size(db), units = "MB"), "\n")
        cat("====================================\n\n")
      }

      # 清理旧文件（可选）
      if(length(existing_files) > 0) {
        if(verbose) cat("清理旧缓存文件...\n")
        for(old_file in existing_files) {
          if(old_file != cache_file) {
            file.remove(old_file)
            if(verbose) cat("  已删除:", basename(old_file), "\n")
          }
        }
        if(verbose) cat("\n")
      }

      return(db)

    }, error = function(e) {
      warning("下载KEGG数据库失败: ", e$message, "\n")

      # 如果下载失败但有旧缓存，使用旧缓存
      if(!is.null(db_cached)) {
        if(verbose) cat("下载失败，使用旧缓存数据\n\n")
        return(db_cached)
      } else {
        stop("无法获取KEGG数据库：下载失败且无本地缓存")
      }
    })
  }
}


#' @title Another matching mode for retrieving KEGG metabolite information based on input metabolite IDs
#' @description
#' This function matches the input metabolite ID to common names in the KEGG database by formatting it,
#' and finally outputs its storage ID in the KEGG database.
#' @param id The metabolite ID for which KEGG information is to be retrieved.
#' @return A data frame with KEGG IDs and matched metabolite names for input metabolite IDs.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' library(dplyr)
#' tax= ps.ms %>% vegan_tax() %>%as.data.frame()
#' id = tax$Metabolite
#' tax2 = ann.kegg2(id)
#' head(tax2)
#' @export

ann.kegg2 <- function(id,
                      cache_dir = "./kegg_cache",
                      max_age_days = 90,
                      force_update = FALSE,
                      verbose = TRUE) {

  # 获取KEGG数据库（带缓存）
  mk <- get_kegg_database(
    cache_dir = cache_dir,
    max_age_days = max_age_days,
    force_update = force_update,
    verbose = verbose
  )

  # 数据预处理
  mk$allMetabolites <- str_to_lower(mk$allMetabolites)

  # 清理输入ID
  id2 <- gsub("[ ][0-9]$", "", id)
  id2 <- str_to_lower(id2)

  # 初始化结果向量
  A <- rep("", length(id))
  B <- rep("", length(id))

  # 准备compound矩阵
  compoundM <- mk[, 1:2]
  colnames(compoundM) <- c("KEGG compound", "common names")
  rownames(compoundM) <- NULL

  # 匹配循环
  if(verbose) cat("开始匹配", length(id), "个代谢物...\n")

  for(i in 1:length(id)) {
    match <- id2[i]

    if(nchar(match) >= 1) {
      target_column <- compoundM[, 2]
      matchM <- match_KEGG(match, target_column, compoundM)

      if(!is.null(matchM[[1]])) {
        A[i] <- matchM[[1]][1, 1]
        B[i] <- matchM[[1]][1, 2]
      }
    }
  }

  # 构建结果
  tax <- data.frame(
    ID = id,
    keggID = A,
    matchname = B,
    stringsAsFactors = FALSE
  )

  # 统计
  matched <- sum(A != "")
  if(verbose) {
    cat("\n========== 匹配统计 ==========\n")
    cat("总数:", length(id), "\n")
    cat("匹配成功:", matched, "\n")
    cat("匹配率:", round(matched/length(id)*100, 2), "%\n")
  }

  return(tax)
}


# 辅助函数：清理KEGG缓存
clear_kegg_cache <- function(cache_dir = "./kegg_cache") {
  if(dir.exists(cache_dir)) {
    unlink(cache_dir, recursive = TRUE)
    cat("缓存已清理:", cache_dir, "\n")
  } else {
    cat("缓存目录不存在\n")
  }
}


# 辅助函数：查看缓存信息
check_kegg_cache <- function(cache_dir = "./kegg_cache") {

  cat("========== KEGG缓存信息 ==========\n")
  cat("缓存目录:", cache_dir, "\n\n")

  # 查找所有缓存文件
  existing_files <- list.files(cache_dir, pattern = "^kegg_database_\\d{8}\\.rds$", full.names = TRUE)

  if(length(existing_files) == 0) {
    cat("缓存文件: 不存在\n")
    cat("状态: 需要下载\n")
  } else {
    cat("找到", length(existing_files), "个缓存文件:\n\n")

    # 按日期排序
    file_dates <- gsub(".*_(\\d{8})\\.rds$", "\\1", basename(existing_files))
    sorted_idx <- order(as.Date(file_dates, format = "%Y%m%d"), decreasing = TRUE)

    for(i in sorted_idx) {
      file <- existing_files[i]
      file_date <- file_dates[i]

      cat("【文件", which(sorted_idx == i), "】\n")
      cat("文件名:", basename(file), "\n")
      cat("日期:", format(as.Date(file_date, format = "%Y%m%d"), "%Y-%m-%d"), "\n")

      # 文件信息
      file_info <- file.info(file)
      cat("大小:", format(file_info$size/1024/1024, digits = 2), "MB\n")

      # 计算已使用天数
      days_old <- as.numeric(Sys.Date() - as.Date(file_date, format = "%Y%m%d"))
      cat("已使用:", days_old, "天\n")

      if(days_old > 90) {
        cat("状态: 已过期\n")
      } else {
        cat("状态: 有效（还有", 90 - days_old, "天过期）\n")
      }

      # 尝试读取记录数
      tryCatch({
        db <- readRDS(file)
        cat("记录数:", nrow(db), "\n")

        # 如果是最新的，标记
        if(i == sorted_idx[1]) {
          cat(">>> 当前使用的版本 <<<\n")
        }
      }, error = function(e) {
        cat("无法读取文件\n")
      })

      cat("\n")
    }
  }

}
