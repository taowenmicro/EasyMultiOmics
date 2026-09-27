# install_EasyMultiOmics.R
# 安装 EasyMultiOmics 及其所有依赖（CRAN、Bioconductor、GitHub）

install.Rtools = function(){

  # 检查并安装 Rtools（Windows 专用）
  if (.Platform$OS.type == "windows") {
    cat("检查 Rtools 安装状态...\n")
    r_version <- paste0(R.version$major, ".", sub("\\..*", "", R.version$minor))
    rtools_needed <- FALSE
    tryCatch({
      if (!("tools" %in% rownames(installed.packages()))) {
        rtools_needed <- TRUE
      } else if (!Sys.which("make") != "") {
        rtools_needed <- TRUE
      }
    }, error = function(e) {
      rtools_needed <- TRUE
    })

    if (rtools_needed) {
      cat("检测到 Rtools 未安装或配置不正确，尝试安装...\n")
      tryCatch({
        # 根据 R 版本选择 Rtools 版本（R >= 4.3.0 使用 Rtools43）
        rtools_url <- if (as.numeric(r_version) >= 4.2) {
          "https://cran.r-project.org/bin/windows/Rtools/rtools45/files/rtools45-6608-6492.exe"
        } else {
          "hhttps://cran.r-project.org/bin/windows/Rtools/rtools40-x86_64.exe"
        }
        temp_file <- tempfile(fileext = ".exe")
        options(timeout = max(6000, getOption("timeout")))
        cat("正在下载 Rtools 从", rtools_url, "...\n")
        download.file(rtools_url, temp_file, mode = "wb")
        cat("正在运行下载的 Rtools 安装程序（需要手动点击安装）：", temp_file, "\n")

        shell(temp_file) # 尝试运行安装程序（可能需管理员权限）
        stop("Rtools 安装程序已下载，请完成安装后重试")
      }, error = function(e) {
        cat("Rtools 下载或安装失败：", e$message, "\n")
        cat("请手动下载并安装 Rtools：", rtools_url, "\n")
        stop("Rtools 安装成功")
      })
    } else {
      cat("Rtools 已安装\n")
    }
  }
}
library(devtools)
install.Rtools()

install.packages("devtools")
install_EasyMultiOmics.liux <- function() {
  # --- 1. 强制设定公共安装路径 ---
  # 这是 Linux 标准的 R 公共库位置，所有用户默认都能读
  public_lib <- "/usr/local/lib/R/site-library"

  cat("==========================================================\n")
  cat("   EasyMultiOmics 多用户部署程序 (强制公共库版) \n")
  cat("==========================================================\n")
  cat("-> 目标公共路径:", public_lib, "\n")

  # 检查并创建目录
  if (!dir.exists(public_lib)) {
    cat("-> 目录不存在，正在创建...\n")
    dir.create(public_lib, recursive = TRUE, mode = "0755")
  }

  # 检查是否有写入权限
  if (file.access(public_lib, mode = 2) != 0) {
    stop("❌ 错误：无法写入目标路径！\n   请务必在 Linux 终端使用 'sudo R' 运行此脚本！")
  }

  # 将公共路径设为首选，确保依赖包也装在这里
  .libPaths(c(public_lib, .libPaths()))

  # 设置超时防止下载中断
  options(timeout = max(1200, getOption("timeout")))

  # --- 2. 基础工具安装 (指定 lib) ---
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cran.r-project.org", lib = public_lib)
  }
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cran.r-project.org", lib = public_lib)
  }
  if (!requireNamespace("pacman", quietly = TRUE)) {
    install.packages("pacman", repos = "https://cran.r-project.org", lib = public_lib)
  }

  # --- 3. 定义依赖列表 ---
  cran_packages <- c(
    "poweRlaw", "randomForest", "fs", "aplot", "Boruta", "broom", "caret",
    "circlize", "classInt", "cluster", "cowplot", "data.tree", "dplyr", "e1071",
    "exactRankTests", "forcats", "ggalluvial", "ggforce", "ggh4x", "ggnewscale",
    "ggplot2", "ggpubr", "ggraph", "ggrepel", "ggsci", "ggstance", "ggtern",
    "ggtree", "ggtreeExtra", "ggupset", "ggVennDiagram", "glmnet", "grid",
    "gridExtra", "gtable", "Hmisc", "igraph", "impute", "MASS", "mboost",
    "MetBrewer", "networkD3", "nlme", "nnet", "packcircles", "patchwork",
    "pheatmap", "plotly", "plspm", "plyr", "ppcor", "preprocessCore", "psych",
    "purrr", "R.utils", "RColorBrewer", "reshape2", "ROCR", "Rtsne", "scales",
    "semPlot", "stats", "stringr", "sunburstR", "tibble", "tidyfst", "tidyr",
    "tidytree", "tidyverse", "treemap", "utils", "vegan", "VennDiagram", "WGCNA",
    "ape", "Ckmeans.1d.dp", "ipred", "picante", "ropls", "rpart", "xgboost"
  )

  bioc_packages <- c(
    "impute","ropls","preprocessCore","ALDEx2", "ANCOMBC", "clusterProfiler",
    "corncob", "DESeq2", "edgeR", "GSVA", "KEGGREST", "lavaan", "limma", "mia",
    "microbiome", "MicrobiotaProcess", "mixOmics", "phyloseq", "sva", "GO.db"
  )

  # --- 4. 安装 CRAN & Bioconductor 包 ---
  # 使用 pacman::p_load 的逻辑手动实现，强制指定 lib

  # 获取已安装包列表 (只看公共库，避免误判)
  installed_pkgs <- installed.packages(lib.loc = public_lib)[, "Package"]

  # 合并需要安装的列表
  all_std <- unique(c(cran_packages, bioc_packages))
  missing_pkgs <- all_std[!all_std %in% installed_pkgs]

  if (length(missing_pkgs) > 0) {
    cat(sprintf("\n[阶段 1] 正在安装 %d 个基础依赖到公共库...\n", length(missing_pkgs)))
    # 这里的关键是 lib = public_lib
    BiocManager::install(missing_pkgs, lib = public_lib, update = FALSE, ask = FALSE)
  } else {
    cat("\n[阶段 1] CRAN/Bioc 基础依赖已就绪。\n")
  }

  # --- 5. 安装 GitHub 包 (强制指定 lib) ---
  github_packages <- c(
    "YuLab-SMU/ggtreeExtra", # 优先修复
    "taowenmicro/ggClusterNet",
    "taowenmicro/EasyMultiOmics.db",
    "taowenmicro/EasyStat",
    "barakbri/dacomp"
  )

  cat("\n[阶段 2] 检查 GitHub 依赖...\n")
  for (pkg in github_packages) {
    pkg_name <- gsub(".*/", "", pkg)
    if (!pkg_name %in% installed.packages(lib.loc = public_lib)[, "Package"]) {
      cat(sprintf("   -> 正在安装: %s ... \n", pkg))
      tryCatch({
        devtools::install_github(pkg, lib = public_lib, upgrade = "never", quiet = TRUE)
        cat("      成功\n")
      }, error = function(e) {
        cat("      失败: ", e$message, "\n")
      })
    }
  }

  # --- 6. 安装主包 EasyMultiOmics ---
  cat("\n[阶段 3] 安装 EasyMultiOmics 主程序...\n")
  devtools::install_github("taowenmicro/EasyMultiOmics",
                           lib = public_lib,
                           upgrade = "never",
                           force = TRUE)

  # --- 7. 关键步骤：修复权限 ---
  cat("\n[阶段 4] 正在修复权限，确保所有用户可用...\n")
  # 这一步调用 Linux 系统命令，把文件夹权限设为 755 (所有人可读可执行)
  system(paste0("chmod -R a+rX ", public_lib))

  # --- 8. 验证 ---
  if (requireNamespace("EasyMultiOmics", lib.loc = public_lib, quietly = TRUE)) {
    cat("\n✅✅✅ 安装成功！ ✅✅✅\n")
    cat("包已安装至公共路径：", public_lib, "\n")
    cat("所有用户现在都可以通过 library(EasyMultiOmics) 使用它了。\n")
  } else {
    stop("❌ 安装似乎失败了，无法从公共路径加载包。")
  }
}
install_EasyMultiOmics <- function() {
  cat("欢迎使用 EasyMultiOmics！开始检查并安装依赖...\n")

  # 设置用户级库路径，避免权限问题
  user_lib <- Sys.getenv("R_LIBS_USER")
  if (!dir.exists(user_lib)) {
    dir.create(user_lib, recursive = TRUE)
    cat("创建用户级库路径：", user_lib, "\n")
  }
  .libPaths(c(user_lib, .libPaths())) # 优先使用用户级路径
  cat("使用库路径：", paste(.libPaths(), collapse = ", "), "\n")

  # 检查路径写入权限
  if (!file.access(user_lib, mode = 2) == 0) {
    cat("警告：用户级库路径", user_lib, "无写入权限！\n")
    stop("请确保", user_lib, "有写入权限，或以管理员身份运行 R")
  }

  # CRAN 依赖
  cran_packages <- c(
    "aplot", "Boruta", "broom", "caret", "circlize", "classInt", "cluster",
    "cowplot", "data.tree", "dplyr", "e1071", "exactRankTests", "forcats",
    "fs", "ggalluvial", "ggforce", "ggh4x", "ggnewscale", "ggplot2", "ggpubr",
    "ggraph", "ggrepel", "ggsci", "ggstance", "ggtern", "ggtree", "ggtreeExtra",
    "ggupset", "ggVennDiagram", "glmnet", "grid", "gridExtra", "gtable", "Hmisc",
    "igraph", "impute", "MASS", "mboost", "MetBrewer", "networkD3", "nlme",
    "nnet", "packcircles", "patchwork", "pheatmap", "plotly", "plspm", "plyr",
    "ppcor", "preprocessCore", "psych", "purrr", "R.utils", "randomForest",
    "RColorBrewer", "reshape2", "ROCR", "Rtsne", "scales", "semPlot", "stats",
    "stringr", "sunburstR", "tibble", "tidyfst", "tidyr", "tidytree",
    "tidyverse", "treemap", "utils", "vegan", "VennDiagram", "WGCNA",
    "ape", "Ckmeans.1d.dp", "ipred", "pacman", "picante", "ropls", "rpart",
    "xgboost"
  )
  tryCatch(if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("正在安装 BiocManager...\n")
    install.packages("BiocManager", repos = "https://cran.r-project.org", quiet = TRUE)
  })
  cran_missing <- cran_packages[!cran_packages %in% installed.packages()[,"Package"]]
  if (length(cran_missing) > 0) {
    cat("检测到缺少 CRAN 依赖：", paste(cran_missing, collapse = ", "), "\n")
    tryCatch({
      cat("正在安装 CRAN 依赖...\n")
      BiocManager::install(cran_missing,  quiet = TRUE)
      cat("CRAN 依赖安装完成！\n")
    }, error = function(e) {
      cat("CRAN 依赖安装失败：", e$message, "\n")
      cat("请手动运行：install.packages(c(", paste(shQuote(cran_missing), collapse = ", "), "))\n")
      stop("依赖安装失败，请检查网络或手动安装")
    })
  }

  # Bioconductor 依赖
  bioc_packages <- c("impute","ropls","preprocessCore",
    "ALDEx2", "ANCOMBC", "clusterProfiler", "corncob",  "DESeq2",
    "edgeR", "GSVA", "KEGGREST", "lavaan", "limma", "mia", "microbiome",
    "MicrobiotaProcess", "mixOmics", "phyloseq", "sva"
  )
  bioc_missing <- bioc_packages[!bioc_packages %in% installed.packages()[,"Package"]]
  if (length(bioc_missing) > 0) {
    cat("检测到缺少 Bioconductor 依赖：", paste(bioc_missing, collapse = ", "), "\n")
    tryCatch({
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        cat("正在安装 BiocManager...\n")
        install.packages("BiocManager", repos = "https://cran.r-project.org", quiet = TRUE)
      }
      cat("正在安装 Bioconductor 依赖...\n")
      BiocManager::install(bioc_missing, ask = FALSE, quiet = TRUE)
      cat("Bioconductor 依赖安装完成！\n")
    }, error = function(e) {
      cat("Bioconductor 依赖安装失败：", e$message, "\n")
      cat("请手动运行：BiocManager::install(c(", paste(shQuote(bioc_missing), collapse = ", "), "))\n")
      stop("依赖安装失败，请检查网络或手动安装")
    })
  }



  tryCatch(if (!requireNamespace("GO.db", quietly = TRUE)) {
    cat("正在安装 BiocManager...\n")
    BiocManager::install("BiocManager", repos = "https://cran.r-project.org", quiet = TRUE)
  }else{
    library("GO.db")
  }
  )

  # GitHub 依赖
  github_packages <- c(
    "taowenmicro/ggClusterNet",
    "taowenmicro/EasyMultiOmics.db",
    "taowenmicro/EasyStat",
    "barakbri/dacomp"
  )
  github_missing <- github_packages[!gsub(".*/", "", github_packages) %in%
                                      installed.packages()[,"Package"]]
  if (length(github_missing) > 0) {
    cat("检测到缺少 GitHub 依赖：", paste(github_missing, collapse = ", "), "\n")
    tryCatch({
      if (!requireNamespace("devtools", quietly = TRUE)) {
        cat("正在安装 devtools...\n")
        BiocManager::install("devtools",  quiet = TRUE)
      }
      cat("正在安装 GitHub 依赖...\n")
      devtools::install_github(github_missing, quiet = TRUE)
      cat("GitHub 依赖安装完成！\n")
    }, error = function(e) {
      cat("GitHub 依赖安装失败：", e$message, "\n")
      cat("请手动运行：devtools::install_github(c(", paste(shQuote(github_missing), collapse = ", "), "))\n")
      stop("依赖安装失败，请检查网络或手动安装")
    })
  }


  # 检查并安装 EasyMultiOmics
  if (length(cran_missing) == 0 && length(bioc_missing) == 0 &&
      length(github_missing) == 0) {
    if (requireNamespace("EasyMultiOmics", quietly = TRUE)) {
      cat("EasyMultiOmics 已安装，准备就绪！\n")
    } else {
      cat("尝试安装 EasyMultiOmics到", user_lib, "...\n")
      tryCatch({
        if (!requireNamespace("devtools", quietly = TRUE)) {
          cat("正在安装 devtools到", user_lib, "...\n")
          install.packages("devtools", repos = "https://cran.r-project.org",
                           lib = user_lib, quiet = TRUE)
        }
        cat("正在从 GitHub 安装 EasyMultiOmics...\n")
        devtools::install_github("taowenmicro/EasyMultiOmics", lib = user_lib,
                                 quiet = TRUE, force = TRUE)
        cat("EasyMultiOmics 安装完成！\n")
      }, error = function(e) {
        cat("EasyMultiOmics 安装失败：", e$message, "\n")
        cat("请手动运行：devtools::install_github('taowenmicro/EasyMultiOmics', lib = '",
            user_lib, "')\n")
        stop("EasyMultiOmics 安装失败，请检查网络或手动安装")
      })
    }
  } else {
    cat("部分依赖安装失败，请检查错误信息或手动安装！\n")
    stop("依赖安装不完整，无法安装 EasyMultiOmics")
  }
}


# 执行安装
install_EasyMultiOmics()


# library("GO.db")

library(EasyMultiOmics)
library(qs)

library(devtools)
install_github("qsbase/qs2")
