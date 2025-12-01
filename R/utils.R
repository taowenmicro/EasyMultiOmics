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
                       width       = NULL,
                       height      = NULL,
                       dpi         = 300,
                       base_width  = 8,
                       base_height = 6,
                       x_per_inch      = 6,
                       facets_per_row  = 4) {

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # estimate number of x-axis categories
  n_x <- NA_integer_
  if (!is.null(p$mapping$x)) {
    xvar <- rlang::as_name(p$mapping$x)
    if (!is.null(p$data) && xvar %in% names(p$data)) {
      n_x <- dplyr::n_distinct(p$data[[xvar]])
    }
  }

  # estimate number of facet panels
  n_facets <- 1L
  gb <- tryCatch(ggplot2::ggplot_build(p),
                 error = function(e) NULL)
  if (!is.null(gb) && !is.null(gb$layout$layout$PANEL)) {
    n_facets <- length(unique(gb$layout$layout$PANEL))
  }

  # width
  if (is.null(width)) {
    add_w <- if (!is.na(n_x)) max(0, n_x / x_per_inch) else 0
    width_final <- base_width + add_w
  } else {
    width_final <- width
  }

  # height
  if (is.null(height)) {
    facet_rows <- ceiling(n_facets / facets_per_row)
    add_h      <- max(0, (facet_rows - 1) * 2)  # +2 inch per extra facet row
    height_final <- base_height + add_h
  } else {
    height_final <- height
  }

  # save PNG
  ggplot2::ggsave(
    filename  = file.path(out_dir, paste0(prefix, ".png")),
    plot      = p,
    width     = width_final,
    height    = height_final,
    dpi       = dpi,
    units     = "in",
    limitsize = FALSE
  )

  # save PDF
  ggplot2::ggsave(
    filename  = file.path(out_dir, paste0(prefix, ".pdf")),
    plot      = p,
    width     = width_final,
    height    = height_final,
    units     = "in"
  )

  invisible(NULL)
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
write_sheet2 <- function(wb, sheet_name, df, row_names = TRUE) {
  stopifnot(inherits(wb, "Workbook"))

  if (sheet_name %in% names(wb)) {
    openxlsx::removeWorksheet(wb, sheet_name)
  }
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, df, rowNames = row_names)

  invisible(NULL)
}
