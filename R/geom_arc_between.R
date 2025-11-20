# ---- geom_arc_between: draw an adjustable arc between two points ----
# Dependencies: ggplot2
#' Draw a curved line (circular arc) between two points
#'
#' @description
#' `geom_arc_between()` draws a circular arc between `(x, y)` and `(xend, yend)`.
#' The arc curvature is controlled by `curv` in [0, 1], where `1` is a semicircle,
#' `0` is a straight line. Use `direction = 1` or `-1` to flip the side of the arc.
#'
#' @param mapping,data,... Standard ggplot2 mapping/data/params passed to `geom_path`.
#'        Aesthetics required: `x`, `y`, `xend`, `yend`.
#' @param curv Numeric in [0, 1]. 0 = straight; 1 = semicircle. Default 1.
#' @param direction Integer, either `1` or `-1`, chooses the side of the arc. Default 1.
#' @param n Integer, number of points used to draw the arc. Default 100.
#' @param lineend Passed to `geom_path`. Default "round".
#' @param na.rm Logical; default FALSE.
#' @param show.legend Logical; default NA.
#' @param inherit.aes Logical; default TRUE.
#'
#' @return A ggplot2 layer (path).
#' @examples
#' library(ggplot2)
#' df <- data.frame(x = 0, y = 0, xend = 4, yend = 0)
#'
#' # semicircle above
#' ggplot(df) +
#'   geom_point(aes(x, y)) +
#'   geom_point(aes(xend, yend), color = "red") +
#'   geom_arc_between(aes(x = x, y = y, xend = xend, yend = yend),
#'                    curv = 1, direction = 1, linewidth = 1) +
#'   coord_equal() + theme_minimal()
#'
#' # gentler arc below
#' ggplot(df) +
#'   geom_arc_between(aes(x, y, xend = xend, yend = yend),
#'                    curv = 0.4, direction = -1, linewidth = 1) +
#'   coord_equal() + theme_minimal()
#'
#' # multiple pairs
#' df2 <- data.frame(
#'   x = c(0, 1), y = c(0, 0),
#'   xend = c(4, 3), yend = c(0, 1),
#'   curv = c(1, 0.6), direction = c(1, -1)
#' )
#' ggplot(df2) +
#'   geom_arc_between(aes(x, y, xend = xend, yend = yend,
#'                        curv = curv, direction = direction),
#'                    linewidth = 1) +
#'   coord_equal() + theme_minimal()
#'
#' @export
geom_arc_between <- function(mapping = NULL, data = NULL, ...,
                             curv = 1,
                             direction = 1,
                             n = 100,
                             lineend = "round",
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE) {
  ggplot2::layer(
    stat = StatArcBetween, data = data, mapping = mapping, geom = ggplot2::GeomPath,
    position = "identity", show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(curv = curv, direction = direction, n = n,
                  lineend = lineend, na.rm = na.rm, ...)
  )
}

#' @keywords internal
StatArcBetween <- ggplot2::ggproto(
  "StatArcBetween", ggplot2::Stat,
  required_aes = c("x", "y", "xend", "yend"),
  compute_panel = function(data, scales,
                           curv = 1,
                           direction = 1,
                           n = 100,
                           lineend = "round",
                           na.rm = FALSE) {
    angle_seq <- function(a1, a2, len, dir_sign) {
      d <- atan2(sin(a2 - a1), cos(a2 - a1))
      if (dir_sign > 0 && d < 0) a2 <- a2 + 2*pi
      if (dir_sign < 0 && d > 0) a2 <- a2 - 2*pi
      seq(a1, a2, length.out = len)
    }

    out_list <- lapply(seq_len(nrow(data)), function(i) {
      x1 <- data$x[i]; y1 <- data$y[i]
      x2 <- data$xend[i]; y2 <- data$yend[i]

      # safer: check columns
      curv_i <- if ("curv" %in% names(data) && !is.na(data$curv[i])) data$curv[i] else curv
      dir_i  <- if ("direction" %in% names(data) && !is.na(data$direction[i])) data$direction[i] else direction

      curv_i <- max(0, min(1, curv_i))
      dir_i  <- ifelse(dir_i >= 0, 1, -1)

      dx <- x2 - x1; dy <- y2 - y1
      L  <- sqrt(dx*dx + dy*dy)
      if (!is.finite(L) || L == 0) return(NULL)

      if (curv_i == 0) {
        return(data.frame(x = c(x1, x2), y = c(y1, y2),
                          PANEL = data$PANEL[i], group = data$group[i]))
      }

      theta <- curv_i * pi
      s_theta2 <- sin(theta/2)
      if (abs(s_theta2) < 1e-6) {
        return(data.frame(x = c(x1, x2), y = c(y1, y2),
                          PANEL = data$PANEL[i], group = data$group[i]))
      }
      R <- L / (2 * s_theta2)
      s <- dir_i * L / (2 * tan(theta/2))

      phi <- atan2(dy, dx)
      ux <- cos(phi); uy <- sin(phi)
      px <- -uy; py <- ux

      mx <- (x1 + x2)/2; my <- (y1 + y2)/2
      cx <- mx + s*px; cy <- my + s*py

      a1 <- atan2(y1 - cy, x1 - cx)
      a2 <- atan2(y2 - cy, x2 - cx)
      ang <- angle_seq(a1, a2, n, dir_i)

      xs <- cx + R*cos(ang)
      ys <- cy + R*sin(ang)

      data.frame(x = xs, y = ys,
                 PANEL = data$PANEL[i], group = data$group[i])
    })

    do.call(rbind, out_list)
  }
)
