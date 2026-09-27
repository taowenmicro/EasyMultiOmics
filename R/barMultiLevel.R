barMultiLevel <- function(ps = NULL, group = "Group", Top = 10,
                                         col.g = NULL, width = 0.6,
                                         label_threshold = 1) {

  library(ggplot2)
  library(ggalluvial)
  library(dplyr)
  library(tidyr)
  library(ggsci)

  # 1. 数据清洗与转化
  tax_mat <- as.matrix(phyloseq::tax_table(ps))
  tax_mat[is.na(tax_mat)] <- "Unknown"
  tax_mat[tax_mat == ""] <- "Unknown"
  phyloseq::tax_table(ps) <- phyloseq::tax_table(tax_mat)

  # 2. 检测并统一列名（兼容标准化前后）
  tax <- as.data.frame(phyloseq::tax_table(ps))
  tax_colnames <- colnames(tax)

  # 检测必需的列是否存在
  required_cols <- c("Kingdom", "Super_Class", "Class", "Sub_Class")
  alt_cols <- c("kingdom", "Super_class", "super_class", "Class", "class",
                "Sub_class", "sub_class")

  # 创建列名映射
  col_mapping <- list(
    Kingdom = c("Kingdom", "kingdom"),
    Super_Class = c("Super_Class", "Super_class", "super_class"),
    Class = c("Class", "class"),
    Sub_Class = c("Sub_Class", "Sub_class", "sub_class")
  )

  # 统一列名
  for (target_name in names(col_mapping)) {
    possible_names <- col_mapping[[target_name]]
    found_name <- NULL

    for (pname in possible_names) {
      if (pname %in% tax_colnames) {
        found_name <- pname
        break
      }
    }

    if (!is.null(found_name) && found_name != target_name) {
      colnames(tax)[colnames(tax) == found_name] <- target_name
      if (exists("verbose") && verbose) {
        message("Renamed column '", found_name, "' to '", target_name, "'")
      }
    }
  }

  # 检查是否所有必需列都存在
  tax_colnames_updated <- colnames(tax)
  missing_cols <- setdiff(required_cols, tax_colnames_updated)

  if (length(missing_cols) > 0) {
    stop("Missing required taxonomy columns: ", paste(missing_cols, collapse = ", "),
         "\nAvailable columns: ", paste(tax_colnames_updated, collapse = ", "),
         "\n\nPlease ensure your data has been annotated with HMDB taxonomy.")
  }

  # 更新 phyloseq 对象
  phyloseq::tax_table(ps) <- phyloseq::tax_table(as.matrix(tax))

  # 3. 转化为相对丰度 (%)
  ps <- ps %>% phyloseq::transform_sample_counts(function(x) {
    x / sum(x, na.rm = TRUE) * 100
  })

  # 4. 提取数据
  otu <- as.data.frame(phyloseq::otu_table(ps))
  if (!phyloseq::taxa_are_rows(ps)) {
    otu <- t(otu)
  }
  otu <- as.data.frame(otu)

  tax <- as.data.frame(phyloseq::tax_table(ps))
  map <- data.frame(phyloseq::sample_data(ps), check.names = FALSE)
  map$Sample <- rownames(map)

  otu$OTU <- rownames(otu)
  tax$OTU <- rownames(tax)

  # 5. 合并长数据
  otu_long <- otu %>%
    tidyr::pivot_longer(cols = -OTU, names_to = "Sample", values_to = "Abundance") %>%
    dplyr::left_join(tax, by = "OTU") %>%
    dplyr::left_join(map[, c("Sample", group)], by = "Sample")

  colnames(otu_long)[colnames(otu_long) == group] <- "Group"

  # 6. 定义层级顺序
  ranks_list <- c("Kingdom", "Super_Class", "Class", "Sub_Class")

  # 7. Top N 筛选逻辑
  get_top_taxa <- function(data, rank, Top) {
    top <- data %>%
      dplyr::group_by(.data[[rank]]) %>%
      dplyr::summarise(Total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(desc(Total)) %>%
      dplyr::filter(.data[[rank]] != "Unknown") %>%
      dplyr::slice(1:Top) %>%
      dplyr::pull(.data[[rank]])
    return(top)
  }

  # 对每一层级进行 Top 标记
  for (rank in ranks_list) {
    top_taxa <- get_top_taxa(otu_long, rank, Top)
    new_col <- paste0(rank, "_plot")
    otu_long[[new_col]] <- ifelse(otu_long[[rank]] %in% top_taxa,
                                  otu_long[[rank]], "Others")
  }

  # 8. 获取所有分组
  groups <- unique(otu_long$Group)

  # 9. 全局颜色分配（所有分组共享）
  all_taxa <- c()
  for (rank in ranks_list) {
    col_name <- paste0(rank, "_plot")
    all_taxa <- c(all_taxa, unique(otu_long[[col_name]]))
  }
  all_taxa <- unique(all_taxa)
  all_taxa <- c(setdiff(all_taxa, "Others"), "Others")

  if (is.null(col.g)) {
    # 为 Kingdom 分配颜色
    kingdom_taxa <- unique(otu_long$Kingdom_plot)
    kingdom_taxa <- kingdom_taxa[kingdom_taxa != "Others"]
    n_kingdom <- length(kingdom_taxa)

    # 为其他层级分配颜色
    super_taxa <- setdiff(unique(otu_long$Super_Class_plot), c("Others", kingdom_taxa))
    class_taxa <- setdiff(unique(otu_long$Class_plot), c("Others", kingdom_taxa, super_taxa))
    sub_taxa <- setdiff(unique(otu_long$Sub_Class_plot), c("Others", kingdom_taxa, super_taxa, class_taxa))

    # 配色方案
    if (n_kingdom > 0) {
      kingdom_colors <- pal_jco()(min(n_kingdom, 10))
      if (n_kingdom > 10) {
        kingdom_colors <- c(kingdom_colors, pal_lancet()(n_kingdom - 10))
      }
      names(kingdom_colors) <- kingdom_taxa
    } else {
      kingdom_colors <- c()
    }

    if (length(super_taxa) > 0) {
      super_colors <- colorRampPalette(c("#E64B35", "#4DBBD5", "#00A087"))(length(super_taxa))
      names(super_colors) <- super_taxa
    } else {
      super_colors <- c()
    }

    if (length(class_taxa) > 0) {
      class_colors <- colorRampPalette(c("#F39B7F", "#8491B4", "#91D1C2"))(length(class_taxa))
      names(class_colors) <- class_taxa
    } else {
      class_colors <- c()
    }

    if (length(sub_taxa) > 0) {
      sub_colors <- colorRampPalette(c("#DC0000", "#7E6148", "#B09C85"))(length(sub_taxa))
      names(sub_colors) <- sub_taxa
    } else {
      sub_colors <- c()
    }

    col.g <- c(kingdom_colors, super_colors, class_colors, sub_colors)
    col.g["Others"] <- "grey80"
    col.g["Unknown"] <- "grey90"
  }

  # 10. 为每个分组创建独立的图
  plot_list <- list()
  flow_data_list <- list()

  for (g in groups) {
    # 筛选当前分组的数据
    group_otu_long <- otu_long %>%
      dplyr::filter(Group == g)

    # 按层级聚合数据
    group_data <- group_otu_long %>%
      dplyr::group_by(Kingdom_plot, Super_Class_plot, Class_plot, Sub_Class_plot) %>%
      dplyr::summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

    # 计算平均丰度
    n_samples <- group_otu_long %>%
      dplyr::pull(Sample) %>%
      unique() %>%
      length()

    group_data$Abundance <- group_data$Abundance / n_samples
    group_data$flow_id <- 1:nrow(group_data)

    # 用于连线的颜色（跟随Kingdom）
    group_data$Kingdom_for_color <- group_data$Kingdom_plot

    # 转换为 ggalluvial 格式
    flow_data <- group_data %>%
      tidyr::pivot_longer(
        cols = c(Kingdom_plot, Super_Class_plot, Class_plot, Sub_Class_plot),
        names_to = "Rank",
        values_to = "Taxon"
      ) %>%
      dplyr::mutate(
        Rank = gsub("_plot", "", Rank),
        Rank = factor(Rank, levels = ranks_list)
      )

    flow_data$Kingdom_color <- factor(flow_data$Kingdom_for_color,
                                      levels = names(col.g)[names(col.g) %in% unique(group_data$Kingdom_plot)])
    flow_data$Taxon <- factor(flow_data$Taxon, levels = names(col.g))

    # 计算 stratum 总丰度
    stratum_abundance <- flow_data %>%
      group_by(Rank, Taxon) %>%
      summarise(Total_Abundance = sum(Abundance), .groups = "drop")

    flow_data <- flow_data %>%
      left_join(stratum_abundance, by = c("Rank", "Taxon"))

    flow_data_list[[g]] <- flow_data

    # 绘制当前分组的桑基图
    p <- ggplot(flow_data, aes(x = Rank, y = Abundance,
                               alluvium = flow_id,
                               stratum = Taxon)) +

      ggalluvial::geom_flow(aes(fill = Kingdom_color), # 连线颜色跟随Kingdom
                            width = width,
                            alpha = 0.6,
                            curve_type = "cubic") +

      ggalluvial::geom_stratum(aes(fill = Taxon), # 节点颜色
                               width = width,
                               color = "white",
                               size = 0.3) +

      geom_text(stat = "stratum",
                aes(label = ifelse(Total_Abundance > label_threshold,
                                   as.character(Taxon), "")),
                size = 3,
                fontface = "bold") +

      scale_fill_manual(
        values = col.g,
        na.value = "grey50"
      ) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_x_discrete(expand = c(0.05, 0.05),
                       labels = c("Kingdom" = "Kingdom",
                                  "Super_Class" = "Super Class",
                                  "Class" = "Class",
                                  "Sub_Class" = "Sub Class")) +
      labs(x = "",
           y = "Mean Relative Abundance (%)",
           title = paste0("Metabolite Hierarchy - ", g)) +  # 标题显示分组名
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, color = "#2C3E50"),
        axis.text.x = element_text(size = 12, face = "bold", color = "black", angle = 0),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10)  # 顶部留更多空间
      )

    plot_list[[g]] <- p
  }

  # 11. 生成汇总统计表
  summary_by_group <- otu_long %>%
    dplyr::select(Sample, Group, Kingdom_plot, Super_Class_plot, Class_plot, Sub_Class_plot, Abundance) %>%
    tidyr::pivot_longer(
      cols = c(Kingdom_plot, Super_Class_plot, Class_plot, Sub_Class_plot),
      names_to = "Rank",
      values_to = "Taxon"
    ) %>%
    dplyr::mutate(Rank = gsub("_plot", "", Rank)) %>%
    dplyr::group_by(Group, Rank, Taxon) %>%
    dplyr::summarise(Mean_Abundance = mean(Abundance), .groups = "drop") %>%
    dplyr::arrange(Group, Rank, desc(Mean_Abundance))

  # 总体统计
  summary_total <- otu_long %>%
    dplyr::select(Sample, Kingdom_plot, Super_Class_plot, Class_plot, Sub_Class_plot, Abundance) %>%
    tidyr::pivot_longer(
      cols = c(Kingdom_plot, Super_Class_plot, Class_plot, Sub_Class_plot),
      names_to = "Rank",
      values_to = "Taxon"
    ) %>%
    dplyr::mutate(Rank = gsub("_plot", "", Rank)) %>%
    dplyr::group_by(Rank, Taxon) %>%
    dplyr::summarise(
      Total_Abundance = sum(Abundance),
      Mean_Abundance = mean(Abundance),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Rank, desc(Total_Abundance))

  # 返回结果
  return(list(
    plot_list = plot_list,        # 多张分组图的列表
    summary = summary_total,       # 总体统计
    summary_by_group = summary_by_group,  # 分组统计
    colors = col.g,                # 颜色映射
    raw_data = otu_long,           # 原始数据
    flow_data_list = flow_data_list  # 每个分组的flow数据
  ))
}
