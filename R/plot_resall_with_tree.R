#' @title Circular Heatmap + Taxonomy Tree from res_all
#' @description
#' Build a taxonomy tree from a phyloseq object and generate a circular
#' consistency heatmap based on `res_all`. Returns both tree and ring plots.
#'
#' @param res_all A list containing:
#'   \describe{
#'     \item{$details}{Per-method results (data.frame), must include columns:
#'       \code{micro}, \code{method}, and optionally \code{adjust.p}}
#'     \item{$summary}{Summary results (data.frame), must include columns:
#'       \code{micro}, \code{n_methods}}
#'   }
#' @param ps A phyloseq object containing taxonomy information.
#' @param topN integer Number of top OTUs/ASVs to display.
#'   If `NULL`, use all in `res_all$summary`.
#' @param type character Layout type, one of `"full"`, `"half"`, or `"3/4"`.
#' @param inner_gap numeric Inner space adjustment (controls the inner radius).
#' @param label_size numeric Font size for labels.
#'
#' @return A list with:
#'   \describe{
#'     \item{ring}{ggplot object, circular heatmap}
#'     \item{tree}{ggtree object, taxonomy tree}
#'     \item{tip_order}{character vector, tip labels order}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' res <- plot_resall_with_tree(res_all, ps, topN = 50, type = "3/4")
#' res$tree   # taxonomy tree
#' res$ring   # ring heatmap
#' res$tip_order
#' }
plot_resall_with_tree <- function(res_all, ps,
                                  topN = NULL,
                                  type = c("full", "half", "3/4"),
                                  inner_gap = 0,
                                  label_size = 3) {
  type <- match.arg(type)

  # --- 1) Extract tables ---
  details_tab <- res_all$details
  summary_tab <- res_all$summary
  res_all$summary$micro %>% unique()
  res_all$summary %>% filter(micro == "s__Terrabacter_sp._MAHUQ.38" )


  if (is.null(topN)) {
    topN <- nrow(summary_tab)
  }

  # --- 2) Add significance flag ---
  details_tab <- details_tab %>%
    mutate(sig = case_when(
      method == "ANCOMII.fast" ~ 1,
      !is.na(adjust.p) & adjust.p < 0.05 ~ 1,
      TRUE ~ 0
    ),
    sig = factor(sig, levels = c(0, 1),
                 labels = c("不显著", "显著")))

  # --- 3) Subset phyloseq ---
  top_OTUs <- summary_tab %>%
    arrange(desc(n_methods)) %>%
    slice_head(n = topN) %>%
    pull(micro)

  ps_sub <- ps %>%
    subset_taxa.wt("OTU", top_OTUs)





  # --- 4) Taxonomy tree ---
  tax <- as.data.frame(phyloseq::tax_table(ps_sub))
  tax$ID <- rownames(tax)

  tax$pathString <- apply(tax, 1, function(x) {
    paste(c("Root", na.omit(x[1:7]), x["ID"]), collapse="/")
  })

  tax_tree <- data.tree::as.Node(tax)
  tax_phylo <- as.phylo(tax_tree)

  # --- 5) Extract tip order from ggtree ---
  p_tree <- ggtree(tax_phylo, layout = "fan", open.angle = 0)
  df <- p_tree$data
  tip_df <- df[df$isTip, ]
  tip_df <- tip_df[order(tip_df$angle), ]
  tip_order <- tip_df$label

  # --- 6) Draw circular heatmap with matching tip order ---
  p_ring <- plot_circular_diffheatmap(details_tab = details_tab,
                                      summary_tab = summary_tab,
                                      topN = topN,
                                      type = type,
                                      inner_gap = inner_gap,
                                      label_size = label_size,
                                      tip_order = NULL)

  return(list(
    ring = p_ring,
    tree = p_tree,
    tip_order = tip_order
  ))
}



#' @title Circular Heatmap for Differential Microbial Analysis
#' @description
#' Generate a circular heatmap that shows the consistency of differential analysis
#' results across multiple methods. Supports full circle, half circle, and 3/4 circle layouts.
#'
#' @param details_tab data.frame A detailed results table with columns:
#'   \describe{
#'     \item{micro}{Microbe/OTU/ASV name}
#'     \item{method}{Name of the differential analysis method}
#'     \item{sig}{Significance marker, can be "显著"/"不显著" or 1/0}
#'   }
#' @param summary_tab data.frame A summary results table with columns:
#'   \describe{
#'     \item{micro}{Microbe/OTU/ASV name (must match details_tab)}
#'     \item{n_methods}{Number of methods that detected this microbe as significant}
#'   }
#' @param topN integer Number of top OTUs/ASVs to display (ranked by n_methods).
#' @param type character Layout type, one of `"full"`, `"half"`, or `"3/4"`.
#' @param inner_gap numeric Inner space adjustment (controls the inner radius of the ring).
#' @param label_size numeric Font size for the outer labels.
#' @param tip_order character Optional vector specifying the OTU/ASV order on the circle.
#'   Typically extracted from a phylogenetic tree tip order. If `NULL`, microbes are
#'   ordered by `n_methods` in summary_tab.
#'
#' @return A ggplot2 object representing the circular heatmap.
#' @export
#'
#' @examples
#' \dontrun{
#' p <- plot_circular_diffheatmap(details_tab, summary_tab, topN = 30, type = "full")
#' print(p)
#' }
plot_circular_diffheatmap <- function(details_tab, summary_tab, topN = 30,
                                      type = c("full", "half", "3/4"),
                                      inner_gap = 0, label_size = 3,
                                      tip_order = NULL) {
  type <- match.arg(type)

  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(ggplot2)
  })
  # summary_tab $micro %>% unique()
  # --- 1) 选 Top OTUs ---
  top_OTUs <- summary_tab %>%
    arrange(desc(n_methods)) %>%
    slice_head(n = topN) %>%
    pull(micro)

  df <- details_tab %>%
    filter(micro %in% top_OTUs) %>%
    mutate(sig = ifelse(sig %in% c("显著", 1), "显著", "不显著"))
  # df$micro%>% unique()
  # --- 2) OTU 顺序 ---
  if (is.null(tip_order)) {
    order_OTUs <- summary_tab %>%
      arrange(desc(n_methods)) %>%
      slice_head(n = topN) %>%
      pull(micro)
  } else {
    order_OTUs <- rev(tip_order)
  }
  df$micro <- factor(df$micro, levels = order_OTUs)

  # --- 3) method → y 轴 ---
  df <- df %>%
    mutate(method_id = as.numeric(factor(method, levels = unique(method))))
  n_methods <- length(unique(df$method_id))

  # --- 4) 根据 type 插 gap ---
  M <- length(order_OTUs)
  gap_n <- switch(
    type,
    "full" = 0L,
    "half" = M,
    "3/4"  = round(M / 3)
  )
  gap_levels <- if (gap_n > 0) paste0("gap_", seq_len(gap_n)) else character(0)

  if (gap_n > 0) {
    df_gap <- expand.grid(
      micro     = gap_levels,
      method_id = unique(df$method_id)
    ) %>%
      mutate(sig = NA_character_)

    df2 <- bind_rows(
      df %>% dplyr::select(micro, method_id, sig),
      df_gap
    )
    df2$micro <- factor(df2$micro, levels = c(order_OTUs, gap_levels))
  } else {
    df2 <- df %>% dplyr::select(micro, method_id, sig)
    df2$micro <- factor(df2$micro, levels = order_OTUs)
  }

  # --- 5) OTU 圆周标签角度与对齐 ---
  df2$micro %>% unique()
  all_levels <- levels(df2$micro)
  n_x        <- length(all_levels)

  label_levels  <- all_levels[!grepl("^gap_", all_levels)]
  idx_label     <- match(label_levels, all_levels)

  angle_rad_lab <- 2 * pi * ((idx_label - 0.5) / n_x)
  angle_deg_lab <- (angle_rad_lab * 180 / pi) %% 360
  text_angle <- ifelse(
    angle_deg_lab > 90 & angle_deg_lab < 270,
    angle_deg_lab + 180, angle_deg_lab
  )
  hjust_val <- ifelse(
    angle_deg_lab > 90 & angle_deg_lab < 270,
    1, 0
  )

  label_df <- data.frame(
    micro      = factor(label_levels, levels = all_levels),
    text_angle = text_angle,
    hjust      = hjust_val,
    stringsAsFactors = FALSE
  )

  # OTU 标签半径
  label_r <- n_methods + 1.2

  # --- 5.1) method (y 轴) 标签：精确放到“图形正上方” ---

  # coord_polar(theta = "x", start = pi/2)
  start_angle <- pi / 2

  idx_all       <- seq_len(n_x)
  # 每个 x 水平对应的极角（度），包含 start 旋转
  angle_rad_all <- 2 * pi * ((idx_all - 0.5) / n_x) + start_angle
  angle_deg_all <- (angle_rad_all * 180 / pi) %% 360

  # 目标角度：图形正上方 = 0° (或 360°)
  target_angle <- 0

  # 圆形差值：保证 359° 和 0° 也能被视为“很近”
  circ_diff <- (angle_deg_all - target_angle + 180) %% 360 - 180
  idx_side  <- which.min(abs(circ_diff))

  side_level <- all_levels[idx_side]

  method_label_df <- df %>%
    dplyr::distinct(method_id, method) %>%
    dplyr::arrange(method_id) %>%
    mutate(
      micro = factor(side_level, levels = all_levels),
      y     = method_id     # 如需整体往里/往外，可以改成 method_id ± 常数
    )

  # --- 6) 绘图 ---
  df2$micro %>% unique()
  p <- ggplot(df2, aes(x = micro, y = method_id, fill = sig)) +
    geom_tile(color = "white", na.rm = TRUE) +
    scale_fill_manual(
      values   = c("不显著" = "grey90", "显著" = "steelblue"),
      na.value = NA,
      name     = "Significance"
    ) +
    coord_polar(theta = "x", start = start_angle, clip = "off") +
    scale_y_continuous(expand = c(0, 0)) +
    ylim(inner_gap, n_methods + 1.5) +
    theme_void(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title      = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.margin     = margin(10, 40, 10, 40),
      axis.text.x     = element_blank(),
      axis.ticks      = element_blank()
    ) +
    labs(
      title = paste0(
        "Circular Heatmap of Differential Analysis (Top ",
        topN, " OTUs)"
      )
    ) +
    # OTU 圆周标签
    geom_text(
      data        = label_df,
      aes(
        x     = micro,
        y     = label_r,
        label = micro,
        angle = -text_angle,
        hjust = hjust
      ),
      inherit.aes = FALSE,
      vjust       = 0,
      size        = label_size
    ) +
    # method (y 轴) 标签：正上方那条径向
    geom_text(
      data        = method_label_df,
      aes(
        x     = as.numeric(micro) + 1,
        y     = y,
        label = method
      ),
      inherit.aes = FALSE,
      size        = label_size,
      angle       = 0,    # 水平文字
      hjust       = 0,    # 左对齐，文字向右伸展
      vjust       = 0.5
    )

  return(p)
}

plot_circular_diffheatmap0 <- function(details_tab, summary_tab, topN = 30,
                                      type = c("full", "half", "3/4"),
                                      inner_gap = 0, label_size = 3,
                                      tip_order = NULL) {
  type <- match.arg(type)

  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(ggplot2)
  })

  # --- 1) Select top OTUs ---
  top_OTUs <- summary_tab %>%
    arrange(desc(n_methods)) %>%
    slice_head(n = topN) %>%
    pull(micro)

  df <- details_tab %>%
    filter(micro %in% top_OTUs) %>%
    mutate(sig = ifelse(sig %in% c("显著", 1), "显著", "不显著"))

  # --- 2) Define order of OTUs ---
  if (is.null(tip_order)) {
    order_OTUs <- summary_tab %>%
      arrange(desc(n_methods)) %>%
      slice_head(n = topN) %>%
      pull(micro)
  } else {
    order_OTUs <- rev(tip_order)
  }
  df$micro <- factor(df$micro, levels = order_OTUs)

  # --- 3) Convert methods to y-axis coordinates ---
  df <- df %>%
    mutate(method_id = as.numeric(factor(method, levels = unique(method))))
  n_methods <- length(unique(df$method_id))

  # --- 4) Insert gaps depending on layout type ---
  M <- length(order_OTUs)
  gap_n <- switch(type,
                  "full" = 0L,
                  "half" = M,
                  "3/4" = round(M / 3))
  gap_levels <- if (gap_n > 0) paste0("gap_", seq_len(gap_n)) else character(0)

  if (gap_n > 0) {
    df_gap <- expand.grid(micro = gap_levels,
                          method_id = unique(df$method_id)) %>%
      mutate(sig = NA_character_)
    df2 <- bind_rows(df %>% dplyr::select(micro, method_id, sig),
                     df_gap)
    df2$micro <- factor(df2$micro, levels = c(order_OTUs, gap_levels))
  } else {
    df2 <- df %>% dplyr::select(micro, method_id, sig)
    df2$micro <- factor(df2$micro, levels = order_OTUs)
  }

  # --- 5) Compute label angles and alignment ---
  all_levels <- levels(df2$micro)
  n_x <- length(all_levels)

  label_levels <- all_levels[!grepl("^gap_", all_levels)]
  idx <- match(label_levels, all_levels)

  angle_rad <- 2 * pi * ((idx - 0.5) / n_x)
  angle_deg <- (angle_rad * 180 / pi) %% 360

  text_angle <- ifelse(angle_deg > 90 & angle_deg < 270,
                       angle_deg + 180, angle_deg)
  hjust_val <- ifelse(angle_deg > 90 & angle_deg < 270, 1, 0)

  label_df <- data.frame(
    micro = factor(label_levels, levels = all_levels),
    text_angle = text_angle,
    hjust = hjust_val,
    stringsAsFactors = FALSE
  )

  # Radius for placing labels
  label_r <- n_methods + 1.2

  # --- 6) Plot circular heatmap ---
  p <- ggplot(df2, aes(x = micro, y = method_id, fill = sig)) +
    geom_tile(color = "white", na.rm = TRUE) +
    scale_fill_manual(values = c("不显著" = "grey90", "显著" = "steelblue"),
                      na.value = NA, name = "Significance") +
    coord_polar(theta = "x", start = pi/2, clip = "off") +
    scale_y_continuous(expand = c(0, 0)) +
    ylim(inner_gap, n_methods + 1.5) +
    theme_void(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.margin = margin(10, 40, 10, 40),
      axis.text.x = element_blank(),
      axis.ticks = element_blank()
    ) +
    labs(title = paste0("Circular Heatmap of Differential Analysis (Top ", topN, " OTUs)")) +
    geom_text(data = label_df,
              aes(x = micro, y = label_r, label = micro,
                  angle = -text_angle, hjust = hjust),
              inherit.aes = FALSE, vjust = 0.5, size = label_size)

  return(p)
}
