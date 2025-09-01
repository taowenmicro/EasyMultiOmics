#' Nature-style publication theme for ggplot2
#'
#' A clean, minimalist theme following Nature journal's visual style
#' with no grid lines and subtle axes, suitable for scientific publications
#'
#' @param base_size Base font size (default: 12)
#' @param base_family Base font family (default: "")
#' @param line_size Line width for axes and ticks (default: 0.3)
#' @param line_color Color for axes and ticks (default: "black")
#' @param title_size_mult Multiplier for title sizes relative to base_size (default: 1.17)
#' @param title_face Font face for titles (default: "bold")
#' @param legend_position Legend position (default: "right")
#'
#' @return A ggplot2 theme object
#' @export
#'
#' @examples
#' library(ggplot2)
#'
#' # Basic usage
#' ggplot(mtcars, aes(wt, mpg, color = factor(cyl))) +
#'   geom_point() +
#'   theme_nature()
#'
#' # With customization
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   theme_nature(base_size = 14, line_size = 0.5)
#'
#' # With legend at bottom
#' ggplot(mtcars, aes(wt, mpg, color = factor(cyl))) +
#'   geom_point() +
#'   theme_nature(legend_position = "bottom")

theme_nature <- function(base_size = 12,
                         base_family = "",
                         line_size = 0.3,
                         line_color = "black",
                         title_size_mult = 1.17,
                         title_face = "bold",
                         legend_position = "right") {

  # Validate inputs
  if (!legend_position %in% c("right", "left", "top", "bottom", "none")) {
    warning("Invalid legend_position. Using 'right'.")
    legend_position <- "right"
  }

  # Calculate relative sizes
  title_size <- base_size * title_size_mult  # ~14 when base is 12

  # Build theme
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
      # Remove all panel elements for clean look
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),

      # Axes styling
      axis.line = ggplot2::element_line(
        linewidth = line_size,  # Changed from 'size' to 'linewidth'
        colour = line_color,
        lineend = "square"
      ),
      axis.ticks = ggplot2::element_line(
        linewidth = line_size,  # Changed from 'size' to 'linewidth'
        colour = line_color
      ),
      axis.ticks.length = ggplot2::unit(0.15, "cm"),

      # Axis titles
      axis.title = ggplot2::element_text(
        size = title_size,
        face = title_face,
        color = line_color
      ),
      axis.title.x = ggplot2::element_text(
        margin = ggplot2::margin(t = 8, unit = "pt")  # Fixed margin syntax
      ),
      axis.title.y = ggplot2::element_text(
        margin = ggplot2::margin(r = 8, unit = "pt")  # Fixed margin syntax
      ),

      # Axis text
      axis.text = ggplot2::element_text(
        size = base_size,
        colour = line_color
      ),

      # Legend styling
      legend.position = legend_position,
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.key.size = ggplot2::unit(1.2, "lines"),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(
        size = base_size,
        color = line_color
      ),
      legend.spacing = ggplot2::unit(0.2, "cm"),

      # Plot title
      plot.title = ggplot2::element_text(
        size = title_size,
        face = title_face,
        hjust = 0.5,
        vjust = 1,
        margin = ggplot2::margin(b = 10, unit = "pt"),  # Fixed margin syntax
        color = line_color
      ),

      # Subtitle and caption
      plot.subtitle = ggplot2::element_text(
        size = base_size,
        hjust = 0.5,
        margin = ggplot2::margin(b = 8, unit = "pt"),  # Fixed margin syntax
        color = line_color
      ),
      plot.caption = ggplot2::element_text(
        size = base_size * 0.8,
        hjust = 1,
        margin = ggplot2::margin(t = 10, unit = "pt"),  # Fixed margin syntax
        color = line_color
      ),

      # Strip styling for facets
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        size = title_size,
        face = title_face,
        color = line_color
      ),

      # Ensure all text uses the base family
      text = ggplot2::element_text(family = base_family),

      # Complete theme
      complete = TRUE
    )
}

#' @rdname theme_nature
#' @export
theme_nature_minimal <- function(base_size = 12,
                                 base_family = "") {
  theme_nature(
    base_size = base_size,
    base_family = base_family,
    line_size = 0.2,
    legend_position = "none"
  )
}

#' @rdname theme_nature
#' @export
theme_nature_presentation <- function(base_size = 16,
                                      base_family = "") {
  theme_nature(
    base_size = base_size,
    base_family = base_family,
    line_size = 0.5,
    title_size_mult = 1.25
  )
}

# Helper function to add common Nature-style modifications
#' Add Nature-style plot modifications
#'
#' @param plot A ggplot object
#' @param expand_x X-axis expansion (default: c(0.02, 0))
#' @param expand_y Y-axis expansion (default: c(0.02, 0))
#'
#' @return Modified ggplot object
#' @export
#'
#' @examples
#' p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
#' nature_style(p)

nature_style <- function(plot,
                         expand_x = c(0.02, 0),
                         expand_y = c(0.02, 0)) {
  plot +
    ggplot2::scale_x_continuous(expand = expand_x) +
    ggplot2::scale_y_continuous(expand = expand_y) +
    theme_nature()
}

# Create the original theme for compatibility
mytheme_nature <- theme_nature()

#' Custom ggplot2 theme for cell-like visualization
#'
#' A clean, minimalist theme suitable for scientific publications
#' with subtle grid lines and professional appearance
#'
#' @param base_size Base font size (default: 12)
#' @param base_family Base font family (default: "")
#' @param grid_color Color for major grid lines (default: "#E5E5E5")
#' @param axis_color Color for axis lines and ticks (default: "#444444")
#' @param text_color_primary Primary text color for titles (default: "#222222")
#' @param text_color_secondary Secondary text color for labels (default: "#333333")
#' @param grid_size Line width for grid lines (default: 0.3)
#' @param axis_size Line width for axes (default: 0.4)
#'
#' @return A ggplot2 theme object
#' @export
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   theme_cell()
#'
#' # With custom colors
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   theme_cell(base_size = 14, grid_color = "#F0F0F0")

theme_cell <- function(base_size = 12,
                       base_family = "",
                       grid_color = "#E5E5E5",
                       axis_color = "#444444",
                       text_color_primary = "#222222",
                       text_color_secondary = "#333333",
                       grid_size = 0.3,
                       axis_size = 0.4) {

  # Calculate relative sizes based on base_size
  size_title <- base_size * 1.33      # ~16 when base is 12
  size_axis_title <- base_size * 1.17 # ~14 when base is 12
  size_axis_text <- base_size         # 12 when base is 12
  size_strip <- base_size * 1.17      # ~14 when base is 12

  # Build theme
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
      # Panel elements
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        color = grid_color,
        linewidth = grid_size  # Changed from 'size' to 'linewidth'
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),

      # Plot background
      plot.background = ggplot2::element_blank(),

      # Axes
      axis.line = ggplot2::element_line(
        linewidth = axis_size,  # Changed from 'size' to 'linewidth'
        colour = axis_color
      ),
      axis.ticks = ggplot2::element_line(
        linewidth = axis_size,  # Changed from 'size' to 'linewidth'
        colour = axis_color
      ),
      axis.title.x = ggplot2::element_text(
        size = size_axis_title,
        face = "bold",
        margin = ggplot2::margin(t = 10, unit = "pt"),  # Fixed margin syntax
        color = text_color_primary
      ),
      axis.title.y = ggplot2::element_text(
        size = size_axis_title,
        face = "bold",
        margin = ggplot2::margin(r = 10, unit = "pt"),  # Fixed margin syntax
        color = text_color_primary
      ),
      axis.text.x = ggplot2::element_text(
        size = size_axis_text,
        colour = text_color_secondary
      ),
      axis.text.y = ggplot2::element_text(
        size = size_axis_text,
        colour = text_color_secondary
      ),

      # Legend
      legend.position = "right",
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(
        size = size_axis_text,
        color = text_color_secondary
      ),
      legend.spacing.x = ggplot2::unit(0.3, "cm"),

      # Plot title
      plot.title = ggplot2::element_text(
        size = size_title,
        face = "bold",
        hjust = 0.5,
        vjust = 1,
        color = text_color_primary
      ),

      # Facet strips
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        size = size_strip,
        face = "bold",
        color = text_color_primary
      ),

      # Ensure all text uses the base family
      text = ggplot2::element_text(family = base_family)
    )
}

# Create a preset variation for dark backgrounds
#' @rdname theme_cell
#' @export
theme_cell_dark <- function(base_size = 12,
                            base_family = "") {
  theme_cell(
    base_size = base_size,
    base_family = base_family,
    grid_color = "#3A3A3A",
    axis_color = "#CCCCCC",
    text_color_primary = "#FFFFFF",
    text_color_secondary = "#E0E0E0"
  ) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "#1E1E1E", color = NA),
      plot.background = ggplot2::element_rect(fill = "#1E1E1E", color = NA)
    )
}

# Usage example with the original style
mytheme_cell <- theme_cell()

# Or with customization
mytheme_cell_custom <- theme_cell(
  base_size = 14,
  grid_color = "#F0F0F0",
  axis_color = "#666666"
)
