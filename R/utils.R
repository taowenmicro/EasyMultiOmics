## utils.R --------------------------------------------------------------
## General helper functions:
## - get_group_cols()     : color mapping for groups
## - save_plot2()         : smart ggplot saver (PNG + PDF)
## - save_circlize_plot() : saver for base / circlize plots (PNG + PDF)
## - write_sheet2()       : unified Excel sheet writer


#' Get a named color vector for groups
#'
#' Generate a named color vector for the supplied groups using
#' palettes from the \pkg{ggsci} package.
#'
#' @param groups Character vector of group names.
#' @param palette Character string, one of \code{"npg"}, \code{"nejm"},
#'   or \code{"lancet"}.
#'
#' @return A named character vector of colors with \code{names(cols) == groups}.
#' @importFrom ggsci pal_npg pal_nejm pal_lancet
#' @export
get_group_cols <- function(groups,
                           palette = c("npg", "nejm", "lancet")) {
  palette <- match.arg(palette)

  groups <- unique(as.character(groups))
  if (length(groups) == 0L) {
    return(character())
  }

  pal_fun <- switch(
    palette,
    npg    = ggsci::pal_npg("nrc"),
    nejm   = ggsci::pal_nejm(),
    lancet = ggsci::pal_lancet()
  )

  cols <- pal_fun(length(groups))
  stats::setNames(cols, groups)
}


#' Save a ggplot object as PNG and PDF with auto-adjusted size
#'
#' This function saves a ggplot object to both PNG and PDF. If \code{width}
#' and/or \code{height} are not specified, they are estimated based on
#' the number of x-axis categories and the number of facet panels.
#'
#' @param p A \code{ggplot} object.
#' @param out_dir Output directory. Will be created if it does not exist.
#' @param prefix File name prefix (without extension).
#' @param width,height Plot width and height in inches. If \code{NULL},
#'   they are estimated automatically.
#' @param dpi Resolution for the PNG file (dots per inch).
#' @param base_width,base_height Base width and height (inches) when estimating size.
#' @param x_per_inch Approximate number of x categories per inch of width.
#' @param facets_per_row Expected number of facet panels per row when estimating height.
#'
#' @return Invisibly returns \code{NULL}; files are written to \code{out_dir}.
#' @export
save_plot2 <- function(p,
                       out_dir,
                       prefix,
                       width        = NULL,
                       height       = NULL,
                       dpi          = 300,
                       base_width   = 8,
                       base_height  = 6,
                       x_per_inch   = 6,
                       facets_per_row = 4,
                       save_png     = TRUE,
                       save_pdf     = TRUE) {

  # 检查对象是否是plot列表
  if (inherits(p, "list")) {
    # 调用save_all_plots函数处理列表中的所有图形
    save_all_plots(plot_list = p, out_dir = out_dir, prefix = prefix,
                   mytheme = NULL, fill_colors = NULL, width = width, height = height)
    return(invisible(NULL))  # 完成后返回，不再继续处理
  }

  # Create directory
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ===================================================================
  # 1. Safely extract x variable name (core fix: supports ~factor(id) etc.)
  # ===================================================================
  n_x <- NA_integer_
  if (!is.null(p$mapping$x)) {
    x_expr <- p$mapping$x

    # Case 1: modern tidy-eval formula like ~factor(id) or ~Species
    if (inherits(x_expr, "formula")) {
      x_expr <- x_expr[[2]]  # strip the ~
    }

    # Now x_expr is either a name/symbol or a call like factor(id)
    x_name <- tryCatch({
      rlang::as_name(x_expr)
    }, error = function(e) NULL)

    # If data exists and variable is present, count distinct levels
    if (!is.null(x_name) && !is.null(p$data) && !is.null(p$data) && x_name %in% names(p$data)) {
      n_x <- dplyr::n_distinct(p$data[[x_name]], na.rm = TRUE)
    }
  }

  # ===================================================================
  # 2. Estimate number of facet panels
  # ===================================================================
  n_facets <- 1L
  gb <- tryCatch(ggplot2::ggplot_build(p), error = function(e) NULL)
  if (!is.null(gb) && !is.null(gb$layout$layout$PANEL)) {
    n_facets <- length(unique(gb$layout$layout$PANEL))
  }

  # ===================================================================
  # 3. Auto-calculate width
  # ===================================================================
  if (is.null(width)) {
    add_w <- if (!is.na(n_x)) max(0, (n_x - 5) / x_per_inch) else 0  # smoother scaling
    width_final <- base_width + add_w

    # Special case: coord_polar() + many categories → needs more space
    if (inherits(p$coordinates, "CoordPolar") && !is.na(n_x) && n_x > 15) {
      width_final <- width_final + 2
    }
  } else {
    width_final <- width
  }

  # ===================================================================
  # 4. Auto-calculate height
  # ===================================================================
  if (is.null(height)) {
    facet_rows <- ceiling(n_facets / facets_per_row)
    add_h <- max(0, (facet_rows - 1) * 2.5)  # 2.5 inches per extra row
    height_final <- base_height + add_h

    # For polar plots, height = width looks better
    if (inherits(p$coordinates, "CoordPolar")) {
      height_final <- width_final
    }
  } else {
    height_final <- height
  }

  # ===================================================================
  # 5. Final save (PNG + PDF)
  # ===================================================================
  if (save_png) {
    ggplot2::ggsave(
      filename  = file.path(out_dir, paste0(prefix, ".png")),
      plot      = p,
      width     = width_final,
      height    = height_final,
      dpi       = dpi,
      units     = "in",
      limitsize = FALSE,
      bg        = "white"
    )
  }

  if (save_pdf) {
    ggplot2::ggsave(
      filename  = file.path(out_dir, paste0(prefix, ".pdf")),
      plot      = p,
      width     = width_final,
      height    = height_final,
      units     = "in",
      bg        = "white"
    )

    message("Saved: ", file.path(out_dir, prefix), ".*   (",
            round(width_final, 1), "\" × ", round(height_final, 1), "\")")

    invisible(NULL)
  }
}


save_plot0 <- function(p,
                       out_dir,
                       prefix,
                       width        = NULL,
                       height       = NULL,
                       dpi          = 300,
                       base_width   = 8,
                       base_height  = 6,
                       x_per_inch   = 6,
                       facets_per_row = 4,
                       save_png     = TRUE,
                       save_pdf     = TRUE) {

  # Create directory
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ===================================================================
  # 1. Safely extract x variable name (core fix: supports ~factor(id) etc.)
  # ===================================================================
  n_x <- NA_integer_
  if (!is.null(p$mapping$x)) {
    x_expr <- p$mapping$x

    # Case 1: modern tidy-eval formula like ~factor(id) or ~Species
    if (inherits(x_expr, "formula")) {
      x_expr <- x_expr[[2]]  # strip the ~
    }

    # Now x_expr is either a name/symbol or a call like factor(id)
    x_name <- tryCatch({
      rlang::as_name(x_expr)
    }, error = function(e) NULL)

    # If data exists and variable is present, count distinct levels
    if (!is.null(x_name) && !is.null(p$data) && !is.null(p$data) && x_name %in% names(p$data)) {
      n_x <- dplyr::n_distinct(p$data[[x_name]], na.rm = TRUE)
    }
  }

  # ===================================================================
  # 2. Estimate number of facet panels
  # ===================================================================
  n_facets <- 1L
  gb <- tryCatch(ggplot2::ggplot_build(p), error = function(e) NULL)
  if (!is.null(gb) && !is.null(gb$layout$layout$PANEL)) {
    n_facets <- length(unique(gb$layout$layout$PANEL))
  }

  # ===================================================================
  # 3. Auto-calculate width
  # ===================================================================
  if (is.null(width)) {
    add_w <- if (!is.na(n_x)) max(0, (n_x - 5) / x_per_inch) else 0  # smoother scaling
    width_final <- base_width + add_w

    # Special case: coord_polar() + many categories → needs more space
    if (inherits(p$coordinates, "CoordPolar") && !is.na(n_x) && n_x > 15) {
      width_final <- width_final + 2
    }
  } else {
    width_final <- width
  }

  # ===================================================================
  # 4. Auto-calculate height
  # ===================================================================
  if (is.null(height)) {
    facet_rows <- ceiling(n_facets / facets_per_row)
    add_h <- max(0, (facet_rows - 1) * 2.5)  # 2.5 inches per extra row
    height_final <- base_height + add_h

    # For polar plots, height = width looks better
    if (inherits(p$coordinates, "CoordPolar")) {
      height_final <- width_final
    }
  } else {
    height_final <- height
  }

  # ===================================================================
  # 5. Final save (PNG + PDF)
  # ===================================================================
  if (save_png) {
    ggplot2::ggsave(
      filename  = file.path(out_dir, paste0(prefix, ".png")),
      plot      = p,
      width     = width_final,
      height    = height_final,
      dpi       = dpi,
      units     = "in",
      limitsize = FALSE,
      bg        = "white"
    )
  }

  if (save_pdf) {
    {
      ggplot2::ggsave(
        filename  = file.path(out_dir, paste0(prefix, ".pdf")),
        plot      = p,
        width     = width_final,
        height    = height_final,
        units     = "in",
        bg        = "white"
      )
    }

    message("Saved: ", file.path(out_dir, prefix), ".*   (",
            round(width_final, 1), "\" × ", round(height_final, 1), "\")")

    invisible(NULL)
  }
}



#' Save a ggplot object as PNG and PDF with auto-adjusted size
#'
#' This function saves a ggplot object to both PNG and PDF. If \code{width}
#' and/or \code{height} are not specified, they are estimated based on
#' the number of x-axis categories and the number of facet panels.
#'
#' @param p A \code{ggplot} object.
#' @param out_dir Output directory. Will be created if it does not exist.
#' @param prefix File name prefix (without extension).
#' @param width,height Plot width and height in inches. If \code{NULL},
#'   they are estimated automatically.
#' @param dpi Resolution for the PNG file (dots per inch).
#' @param base_width,base_height Base width and height (inches) when estimating size.
#' @param x_per_inch Approximate number of x categories per inch of width.
#' @param facets_per_row Expected number of facet panels per row when estimating height.
#'
#' @return Invisibly returns \code{NULL}; files are written to \code{out_dir}.
#' @export
save_plot3 <- function(p,
                       out_dir,
                       prefix,
                       width       = NULL,
                       height      = NULL,
                       dpi         = 300,
                       base_width  = 8,
                       base_height = 6,
                       x_per_inch      = 6,
                       facets_per_row  = 4) {

  # 允许 ggplot / gg / aplot / grob
  if (!inherits(p, c("gg", "ggplot", "aplot", "grob"))) {
    stop("`p` 必须是 ggplot / aplot / grob 对象，当前 class: ",
         paste(class(p), collapse = ", "))
  }

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  ## ---------- 选一个代表性的 ggplot 用来估计宽高 ----------
  p_est <- NULL
  if (inherits(p, c("gg", "ggplot"))) {
    p_est <- p
  } else if (inherits(p, "aplot")) {
    # aplot 本质是一个包含多个 ggplot/grob 的列表
    for (i in seq_along(p)) {
      if (inherits(p[[i]], c("gg", "ggplot"))) {
        p_est <- p[[i]]
        break
      }
    }
  }

  # 默认值
  n_x      <- NA_integer_
  n_facets <- 1L

  # 只有在能找到一个 ggplot 子图的前提下，才去估计
  if (!is.null(p_est)) {

    # 估计 x 轴类别数
    if (!is.null(p_est$mapping) && !is.null(p_est$mapping$x)) {
      xvar <- rlang::as_name(p_est$mapping$x)
      if (!is.null(p_est$data) && xvar %in% names(p_est$data)) {
        n_x <- dplyr::n_distinct(p_est$data[[xvar]])
      }
    }

    # 估计 facet panel 数
    gb <- tryCatch(ggplot2::ggplot_build(p_est),
                   error = function(e) NULL)
    if (!is.null(gb) && !is.null(gb$layout$layout$PANEL)) {
      n_facets <- length(unique(gb$layout$layout$PANEL))
    }
  }

  ## ---------- 计算宽度 ----------
  if (is.null(width)) {
    add_w <- if (!is.na(n_x)) max(0, n_x / x_per_inch) else 0
    width_final <- base_width + add_w
  } else {
    width_final <- width
  }

  ## ---------- 计算高度 ----------
  if (is.null(height)) {
    facet_rows <- ceiling(n_facets / facets_per_row)
    add_h      <- max(0, (facet_rows - 1) * 2)  # 每多一行 facet+2 inch
    height_final <- base_height + add_h
  } else {
    height_final <- height
  }

  ## ---------- 存 PNG ----------
  ggplot2::ggsave(
    filename  = file.path(out_dir, paste0(prefix, ".png")),
    plot      = p,
    width     = width_final,
    height    = height_final,
    dpi       = dpi,
    units     = "in",
    limitsize = FALSE
  )

  ## ---------- 存 PDF ----------
  ggplot2::ggsave(
    filename  = file.path(out_dir, paste0(prefix, ".pdf")),
    plot      = p,
    width     = width_final,
    height    = height_final,
    units     = "in"
  )

  invisible(list(width = width_final, height = height_final))
}


#' Save circlize or base graphics plots as PNG and PDF
#'
#' Evaluate an expression that draws a plot (e.g. using base graphics or
#' \pkg{circlize}) and save it as both PNG and PDF files.
#'
#' @param expr An expression that produces a plot when evaluated.
#' @param out_dir Output directory.
#' @param prefix File name prefix (without extension).
#' @param width,height Device width and height (inches).
#' @param res Resolution for the PNG device (dots per inch).
#'
#' @return Invisibly returns \code{NULL}; files are written to \code{out_dir}.
#' @importFrom grDevices png pdf dev.off
#' @export
save_circlize_plot <- function(expr,
                               out_dir,
                               prefix,
                               width  = 10,
                               height = 10,
                               res    = 300) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # PNG
  grDevices::png(
    filename = file.path(out_dir, paste0(prefix, ".png")),
    width    = width,
    height   = height,
    units    = "in",
    res      = res
  )
  eval(substitute(expr), envir = parent.frame())
  grDevices::dev.off()

  # PDF
  grDevices::pdf(
    file   = file.path(out_dir, paste0(prefix, ".pdf")),
    width  = width,
    height = height
  )
  eval(substitute(expr), envir = parent.frame())
  grDevices::dev.off()

  invisible(NULL)
}


#' Write (and overwrite) a worksheet in an Excel workbook
#'
#' Helper wrapper around \pkg{openxlsx} to always overwrite an existing
#' sheet before writing new data.
#'
#' @param wb An \code{openxlsx::Workbook} object.
#' @param sheet_name Name of the worksheet to create/overwrite.
#' @param df A data frame to write.
#' @param row_names Logical, whether to write row names (default: \code{TRUE}).
#'
#' @return Invisibly returns \code{NULL}; modifies \code{wb} in place.
#' @importFrom openxlsx removeWorksheet addWorksheet writeData
#' @export
write_sheet2 <- function(wb,
                         sheet_name,
                         data,
                         row_names = TRUE,
                         suffix_separator = "_",
                         remove_plot_suffix = TRUE,
                         verbose = FALSE) {

  # 检查输入
  if (!inherits(wb, "Workbook")) {
    stop("'wb' must be an openxlsx Workbook object")
  }

  if (!is.character(sheet_name) || length(sheet_name) != 1) {
    stop("'sheet_name' must be a single character string")
  }

  if (is.data.frame(data)) {
    # 验证数据
    if (is.null(data) || nrow(data) == 0) {
      if (verbose) {
        message("Skipping empty data frame for sheet '", sheet_name, "'")
      }
      return(invisible(NULL))
    }

    # 检查sheet名称长度（Excel限制31字符）
    if (nchar(sheet_name) > 31) {
      warning("Sheet name '", sheet_name, "' exceeds 31 characters. Truncating.")
      sheet_name <- substr(sheet_name, 1, 31)
    }

    # 删除已存在的sheet
    if (sheet_name %in% names(wb)) {
      if (verbose) {
        message("Removing existing sheet '", sheet_name, "'")
      }
      openxlsx::removeWorksheet(wb, sheet_name)
    }

    # 添加新sheet
    openxlsx::addWorksheet(wb, sheet_name)

    # 写入数据
    if (verbose) {
      message("Writing ", nrow(data), " rows to sheet '", sheet_name, "'")
    }

    openxlsx::writeData(
      wb = wb,
      sheet = sheet_name,
      x = data,
      rowNames = row_names
    )

    return(invisible(NULL))
  }

  if (is.list(data)) {
    # 验证输入
    if (is.null(data) || length(data) == 0) {
      if (verbose) {
        message("Empty data list provided, nothing to write")
      }
      return(invisible(NULL))
    }

    # 确保列表有名称
    if (is.null(names(data))) {
      names(data) <- paste0("Sheet", seq_along(data))
      warning("Unnamed list provided, using default names: ",
              paste(names(data), collapse = ", "))
    }

    # 统计信息
    total_sheets <- length(data)
    successful_writes <- 0
    failed_writes <- 0
    skipped_writes <- 0

    if (verbose) {
      message("\n========== Writing ", total_sheets, " sheets ==========")
    }

    # 遍历列表中的每个数据框
    for (i in seq_along(data)) {
      name <- names(data)[i]
      df <- data[[i]]

      # 跳过空数据
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
        if (verbose) {
          message("[", i, "/", total_sheets, "] Skipping empty/invalid data: ", name)
        }
        skipped_writes <- skipped_writes + 1
        next
      }

      # 清理名称
      clean_name <- name
      if (remove_plot_suffix) {
        clean_name <- gsub("\\.plot$", "", clean_name)
        clean_name <- gsub("_plot$", "", clean_name)
      }

      # 构建完整的sheet名称
      full_sheet_name <- paste0(sheet_name, suffix_separator, clean_name)

      # 处理过长的名称
      if (nchar(full_sheet_name) > 31) {
        # 尝试缩短prefix
        max_prefix_len <- 31 - nchar(clean_name) - nchar(suffix_separator) - 1
        if (max_prefix_len > 0) {
          truncated_prefix <- substr(sheet_name, 1, max_prefix_len)
          full_sheet_name <- paste0(truncated_prefix, suffix_separator, clean_name)
        } else {
          # 如果还是太长，直接截断
          full_sheet_name <- substr(full_sheet_name, 1, 31)
        }

        if (verbose) {
          message("  Note: Sheet name truncated to: ", full_sheet_name)
        }
      }

      # 写入sheet
      tryCatch({
        # 删除已存在的sheet
        if (full_sheet_name %in% names(wb)) {
          openxlsx::removeWorksheet(wb, full_sheet_name)
        }

        # 添加新sheet
        openxlsx::addWorksheet(wb, full_sheet_name)

        # 写入数据
        openxlsx::writeData(
          wb = wb,
          sheet = full_sheet_name,
          x = df,
          rowNames = row_names
        )

        if (verbose) {
          message("[", i, "/", total_sheets, "] ✓ ", full_sheet_name,
                  " (", nrow(df), " rows)")
        }

        successful_writes <- successful_writes + 1

      }, error = function(e) {
        warning("Failed to write sheet '", full_sheet_name, "': ", e$message)
        if (verbose) {
          message("[", i, "/", total_sheets, "] ✗ ", full_sheet_name,
                  " - ERROR: ", e$message)
        }
        failed_writes <<- failed_writes + 1
      })
    }



    return(invisible(NULL))
  }

  stop("'data' must be either a data frame or a list of data frames")
}


#' @export
write_sheet0 <- function(wb, sheet_name, df, row_names = TRUE) {
  stopifnot(inherits(wb, "Workbook"))

  if (sheet_name %in% names(wb)) {
    openxlsx::removeWorksheet(wb, sheet_name)
  }
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, df, rowNames = row_names)

  invisible(NULL)
}

#' @export
save_all_plots <- function(plot_list,
                           out_dir,
                           prefix = "plot",
                           mytheme = NULL,
                           fill_colors = NULL,
                           width = 10,
                           height = 8) {

  # 如果列表为空，直接返回
  if (is.null(plot_list) || length(plot_list) == 0) {
    message("The plot list is empty. Nothing to save.")
    return(invisible(NULL))
  }

  # 创建输出目录
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # 获取列表的名字 (比如 "OE.WT", "KO.OE")
  list_names <- names(plot_list)

  for (i in seq_along(plot_list)) {
    # 1. 获取绘图对象
    p <- plot_list[[i]]

    # 2. 应用主题 (如果提供了)
    if (!is.null(mytheme)) {
      p <- p + mytheme
    }

    # 3. 应用颜色 (如果提供了)
    if (!is.null(fill_colors)) {
      p <- p +
        guides(fill = guide_legend(title = NULL)) +
        scale_fill_manual(values = fill_colors)
    }

    # 4. 动态生成文件名
    # 如果列表有名字且不为空，使用名字；否则使用数字索引
    if (!is.null(list_names) && list_names[i] != "") {
      suffix <- list_names[i]
    } else {
      suffix <- as.character(i)
    }

    # 拼接前缀和后缀，例如: "pathway_enrich" + "_" + "OE.WT"
    plot_name <- paste0(prefix, "_", suffix)

    # 5. 保存
    save_plot2(p, out_dir = out_dir, prefix = plot_name)
  }

  invisible(NULL)
}

#' @export
write_all_sheets <- function(wb, data_list, prefix) {
  if (is.null(data_list) || length(data_list) == 0) return(invisible(NULL))

  # 遍历所有数据
  for (name in names(data_list)) {
    data <- data_list[[name]]
    if (is.null(data) || nrow(data) == 0) next

    # 构建sheet名称：prefix_name（去掉.plot后缀）
    sheet_name <- paste0(prefix, "_", gsub("\\.plot$", "", name))

    # 写入sheet
    tryCatch({
      write_sheet2(wb, sheet_name, data)
    }, error = function(e) {
      warning("Failed to write sheet '", sheet_name, "': ", e$message)
    })
  }

  invisible(wb)
}
