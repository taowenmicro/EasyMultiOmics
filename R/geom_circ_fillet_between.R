#' Circular/elliptical filleted elbow between two points (A->P->B, ∠APB=90°)
#'
#' @description
#' Connect (x, y) and (xend, yend) by bending at a pivot P on the circle whose
#' diameter is AB (Thales theorem ⇒ ∠APB=90°). Replace the right-angle at P
#' with a quarter circle/ellipse that is tangent to both straight segments.
#'
#' @param mapping,data,... ggplot2 standard. Required aesthetics: `x`, `y`, `xend`, `yend`.
#' @param t Numeric in [0,1]. Pivot position along the chosen semicircle from A to B.
#' @param side "ccw" or "cw": choose counter-clockwise / clockwise semicircle (viewed from A to B).
#' @param corner_rx,corner_ry Numeric. Fillet radii along PA and PB, respectively.
#'   Set equal for a circular fillet; set unequal to control ellipticity.
#' @param n Integer. Number of points to draw the fillet arc. Default 30.
#' @param na.rm,show.legend,inherit.aes Passed to ggplot2.
#'
#' @return A ggplot2 layer (GeomPath).
#' @examples
#' library(ggplot2)
#' df <- data.frame(x=0, y=0, xend=4, yend=2)
#'
#' # 圆角（r=0.4），逆时针半圆
#' ggplot(df) +
#'   geom_circ_fillet_between(aes(x,y,xend=xend,yend=yend),
#'                            t=0.5, side="ccw", corner_rx=0.4, corner_ry=0.4,
#'                            linewidth=1) +
#'   coord_equal() + theme_minimal()
#'
#' # 椭圆角（rx != ry），顺时针半圆
#' ggplot(df) +
#'   geom_circ_fillet_between(aes(x,y,xend=xend,yend=yend),
#'                            t=0.3, side="cw", corner_rx=0.6, corner_ry=0.2,
#'                            linewidth=1) +
#'   coord_equal() + theme_minimal()
#'
#' @export
geom_circ_fillet_between <- function(mapping = NULL, data = NULL, ...,
                                     t = 0.5, side = c("ccw","cw"),
                                     corner_rx = 0.3, corner_ry = 0.3,
                                     n = 30,
                                     na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  side <- match.arg(side)
  ggplot2::layer(
    stat = StatCircFilletBetween, data = data, mapping = mapping,
    geom = ggplot2::GeomPath, position = "identity",
    show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(t = t, side = side, corner_rx = corner_rx, corner_ry = corner_ry,
                  n = n, na.rm = na.rm, ...)
  )
}

#' @keywords internal
StatCircFilletBetween <- ggplot2::ggproto(
  "StatCircFilletBetween", ggplot2::Stat,
  required_aes = c("x","y","xend","yend"),
  compute_panel = function(data, scales, t=0.5, side="ccw",
                           corner_rx=0.3, corner_ry=0.3, n=30) {

    side_sign <- if (identical(side, "ccw")) +1 else -1

    build_row <- function(i) {
      x1 <- data$x[i]; y1 <- data$y[i]
      x2 <- data$xend[i]; y2 <- data$yend[i]

      grp <- if ("group" %in% names(data)) data$group[i] else i
      style_cols <- intersect(c("colour","alpha","linetype","linewidth"), names(data))
      style_vals <- lapply(style_cols, function(cl) data[[cl]][i]); names(style_vals) <- style_cols

      # Circle with AB as diameter
      mx <- (x1 + x2)/2; my <- (y1 + y2)/2
      R  <- sqrt((x1 - mx)^2 + (y1 - my)^2)
      if (!is.finite(R) || R == 0) {
        out <- data.frame(x=c(x1,x2), y=c(y1,y2), PANEL=data$PANEL[i], group=grp)
        for (cl in style_cols) out[[cl]] <- style_vals[[cl]]
        return(out)
      }

      # Pivot P on chosen semicircle
      phiA <- atan2(y1 - my, x1 - mx)
      tt   <- if ("t" %in% names(data) && !is.na(data$t[i])) data$t[i] else t
      tt   <- max(0, min(1, tt))
      phiP <- phiA + side_sign * tt * pi
      Px <- mx + R * cos(phiP)
      Py <- my + R * sin(phiP)

      # Unit vectors along legs (P->A and P->B), guaranteed orthogonal
      u1x <- (x1 - Px); u1y <- (y1 - Py); L1 <- sqrt(u1x^2 + u1y^2); u1x <- u1x / L1; u1y <- u1y / L1
      u2x <- (x2 - Px); u2y <- (y2 - Py); L2 <- sqrt(u2x^2 + u2y^2); u2x <- u2x / L2; u2y <- u2y / L2

      # Clamp radii so tangency points remain on segments
      rx <- max(0, min(corner_rx, 0.99 * L1))
      ry <- max(0, min(corner_ry, 0.99 * L2))

      if (rx == 0 && ry == 0) {
        # Sharp right angle
        xs <- c(x1, Px, x2); ys <- c(y1, Py, y2)
      } else {
        if (rx == 0 || ry == 0) { # one-sided fillet degenerates to a kink + short arc -> treat as sharp
          xs <- c(x1, Px, x2); ys <- c(y1, Py, y2)
        } else {
          # Tangency points on the legs (back off from P toward A by rx; from P toward B by ry)
          E1x <- Px + u1x * rx; E1y <- Py + u1y * rx
          E2x <- Px + u2x * ry; E2y <- Py + u2y * ry

          # Ellipse center so that arc is tangent to both legs:
          # C = P + u1*rx + u2*ry
          Cx <- Px + u1x * rx + u2x * ry
          Cy <- Py + u1y * rx + u2y * ry

          # Parameterize quarter ellipse from E1 to E2:
          # Q(θ) = C - rx * u1 * sinθ - ry * u2 * cosθ, θ ∈ [0, π/2]
          theta <- seq(0, pi/2, length.out = max(2, n))
          arc_x <- Cx - rx * u1x * sin(theta) - ry * u2x * cos(theta)
          arc_y <- Cy - rx * u1y * sin(theta) - ry * u2y * cos(theta)

          xs <- c(x1, E1x, arc_x, E2x, x2)
          ys <- c(y1, E1y, arc_y, E2y, y2)
        }
      }

      out <- data.frame(x = xs, y = ys, PANEL = data$PANEL[i], group = grp)
      for (cl in style_cols) out[[cl]] <- style_vals[[cl]]
      out
    }

    do.call(rbind, lapply(seq_len(nrow(data)), build_row))
  }
)
