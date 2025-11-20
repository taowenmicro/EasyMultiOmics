#' Rectangle-routed polyline with two filleted right-angle bends (A -> P1 -> P2 -> B)
#'
#' @description
#' Connect (x, y) to (xend, yend) via two orthogonal bends at imaginary points P1 and P2
#' so that A, P1, P2, B form an axis-aligned rectangle. The path is A -> P1 -> P2 -> B.
#' At each bend (P1, P2) the 90° corner is replaced by a quarter circle/ellipse that is
#' tangent to both incident segments.
#'
#' @param mapping,data,... ggplot2 standard arguments. Required aesthetics: `x`, `y`, `xend`, `yend`.
#' @param axis Character, either `"x"` (vertical bus x = constant) or `"y"` (horizontal bus y = constant).
#' @param t Numeric in [0, 1]. Relative position of the bus between A and B along the chosen axis.
#'   Ignored if `bus` is provided. Default `0.5` (middle).
#' @param bus Numeric (optional). Absolute coordinate of the bus: if `axis="x"`, this is `xm`;
#'   if `axis="y"`, this is `ym`. Overrides `t` when not `NULL`.
#' @param rx1,ry1 Numeric. Fillet radii at the first bend (P1), along incoming and outgoing legs.
#'   Set `rx1 = ry1` for circular fillet; unequal values give an elliptical fillet. Default `0.3, 0.3`.
#' @param rx2,ry2 Numeric. Fillet radii at the second bend (P2). Default `0.3, 0.3`.
#' @param n Integer. Number of points to draw each fillet arc. Default `30`.
#' @param na.rm,show.legend,inherit.aes Passed to ggplot2.
#'
#' @return A ggplot2 layer (`GeomPath`) drawing A -> P1 -> P2 -> B with two rounded right angles.
#'
#' @examples
#' library(ggplot2)
#' df <- data.frame(x=0, y=0, xend=6, yend=3)
#'
#' # 垂直中线（axis="x"），中线在 A 与 B 的 x 中点，两个圆角半径相同（圆角）
#' ggplot(df) +
#'   geom_rect_fillet_between(aes(x, y, xend = xend, yend = yend),
#'                            axis = "x", t = 0.5,
#'                            rx1 = 0.5, ry1 = 0.5, rx2 = 0.5, ry2 = 0.5,
#'                            linewidth = 1) +
#'   coord_equal() + theme_minimal()
#'
#' # 水平中线（axis="y"），指定绝对 y=1.2，使用椭圆角（扁/瘦可由 rx != ry 控制）
#' ggplot(df) +
#'   geom_rect_fillet_between(aes(x, y, xend = xend, yend = yend),
#'                            axis = "y", bus = 1.2,
#'                            rx1 = 0.6, ry1 = 0.2, rx2 = 0.4, ry2 = 0.6,
#'                            linewidth = 1) +
#'   coord_equal() + theme_minimal()
#'
#' @export
geom_rect_fillet_between <- function(mapping = NULL, data = NULL, ...,
                                     axis = c("x","y"),
                                     t = 0.5, bus = NULL,
                                     rx1 = 0.3, ry1 = 0.3,
                                     rx2 = 0.3, ry2 = 0.3,
                                     n = 30,
                                     na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  axis <- match.arg(axis)
  ggplot2::layer(
    stat = StatRectFilletBetween, data = data, mapping = mapping,
    geom = ggplot2::GeomPath, position = "identity",
    show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(axis = axis, t = t, bus = bus,
                  rx1 = rx1, ry1 = ry1, rx2 = rx2, ry2 = ry2,
                  n = n, na.rm = na.rm, ...)
  )
}

# internal
StatRectFilletBetween <- ggplot2::ggproto(
  "StatRectFilletBetween", ggplot2::Stat,
  required_aes = c("x","y","xend","yend"),
  compute_panel = function(data, scales,
                           axis = "x", t = 0.5, bus = NULL,
                           rx1 = 0.3, ry1 = 0.3, rx2 = 0.3, ry2 = 0.3,
                           n = 30) {

    # helper: build one quarter ellipse fillet at corner P between Pin->P and P->Pout
    fillet_parts <- function(Pin, P, Pout, rx, ry, n) {
      u1 <- Pin - P; L1 <- sqrt(sum(u1^2)); if (L1 == 0) return(NULL); u1 <- u1 / L1
      u2 <- Pout - P; L2 <- sqrt(sum(u2^2)); if (L2 == 0) return(NULL); u2 <- u2 / L2
      # clamp radii
      rx <- max(0, min(rx, 0.99 * L1))
      ry <- max(0, min(ry, 0.99 * L2))
      if (rx == 0 || ry == 0) {
        # no rounding: return split points
        return(list(Ein = P, arc_x = numeric(0), arc_y = numeric(0), Eout = P))
      }
      Ein <- P + u1 * rx
      Eout <- P + u2 * ry
      C <- P + u1 * rx + u2 * ry
      theta <- seq(0, pi/2, length.out = max(2, n))
      arc <- t(sapply(theta, function(th) C - rx * u1 * sin(th) - ry * u2 * cos(th)))
      list(Ein = Ein, arc_x = arc[,1], arc_y = arc[,2], Eout = Eout)
    }

    out_list <- lapply(seq_len(nrow(data)), function(i) {
      x1 <- data$x[i]; y1 <- data$y[i]
      x2 <- data$xend[i]; y2 <- data$yend[i]

      grp <- if ("group" %in% names(data)) data$group[i] else i
      style_cols <- intersect(c("colour","alpha","linetype","linewidth"), names(data))
      style_vals <- lapply(style_cols, function(cl) data[[cl]][i]); names(style_vals) <- style_cols

      # Degenerate: A==B -> straight segment
      if (!is.finite(x1 + y1 + x2 + y2) || (abs(x1 - x2) + abs(y1 - y2) < 1e-12)) {
        df <- data.frame(x = c(x1, x2), y = c(y1, y2),
                         PANEL = data$PANEL[i], group = grp)
        for (cl in style_cols) df[[cl]] <- style_vals[[cl]]
        return(df)
      }

      # Determine bus coordinate
      if (!is.null(bus)) {
        xm <- bus; ym <- bus
      } else {
        tt <- max(0, min(1, if ("t" %in% names(data) && !is.na(data$t[i])) data$t[i] else t))
        xm <- (1 - tt) * x1 + tt * x2
        ym <- (1 - tt) * y1 + tt * y2
      }

      # Build rectangle route via P1, P2
      if (axis == "x") {
        # vertical bus x = xm
        P1 <- c(xm, y1); P2 <- c(xm, y2)
        # straight parts before/after fillets
        # A -> P1 corner fillet
        f1 <- fillet_parts(Pin = c(x1, y1), P = P1, Pout = P2, rx = rx1, ry = ry1, n = n)
        # middle segment: from f1$Eout to P2 (until start of second fillet)
        # P2 corner fillet
        f2 <- fillet_parts(Pin = P1, P = P2, Pout = c(x2, y2), rx = rx2, ry = ry2, n = n)

        # assemble polyline: A -> (to Ein1) -> arc1 -> (to Ein2) -> arc2 -> (to B)
        # segment A -> Ein1
        seg1_x <- c(x1, f1$Ein[1]); seg1_y <- c(y1, f1$Ein[2])
        # arc1
        arc1_x <- c(f1$Ein[1], f1$arc_x, f1$Eout[1]); arc1_y <- c(f1$Ein[2], f1$arc_y, f1$Eout[2])
        # straight to Ein2
        seg2_x <- c(f1$Eout[1], f2$Ein[1]); seg2_y <- c(f1$Eout[2], f2$Ein[2])
        # arc2
        arc2_x <- c(f2$Ein[1], f2$arc_x, f2$Eout[1]); arc2_y <- c(f2$Ein[2], f2$arc_y, f2$Eout[2])
        # to B
        seg3_x <- c(f2$Eout[1], x2); seg3_y <- c(f2$Eout[2], y2)

        xs <- c(seg1_x, arc1_x[-1], seg2_x[-1], arc2_x[-1], seg3_x[-1])
        ys <- c(seg1_y, arc1_y[-1], seg2_y[-1], arc2_y[-1], seg3_y[-1])

      } else {
        # axis == "y": horizontal bus y = ym
        P1 <- c(x1, ym); P2 <- c(x2, ym)

        f1 <- fillet_parts(Pin = c(x1, y1), P = P1, Pout = P2, rx = rx1, ry = ry1, n = n)
        f2 <- fillet_parts(Pin = P1, P = P2, Pout = c(x2, y2), rx = rx2, ry = ry2, n = n)

        seg1_x <- c(x1, f1$Ein[1]); seg1_y <- c(y1, f1$Ein[2])
        arc1_x <- c(f1$Ein[1], f1$arc_x, f1$Eout[1]); arc1_y <- c(f1$Ein[2], f1$arc_y, f1$Eout[2])
        seg2_x <- c(f1$Eout[1], f2$Ein[1]); seg2_y <- c(f1$Eout[2], f2$Ein[2])
        arc2_x <- c(f2$Ein[1], f2$arc_x, f2$Eout[1]); arc2_y <- c(f2$Ein[2], f2$arc_y, f2$Eout[2])
        seg3_x <- c(f2$Eout[1], x2); seg3_y <- c(f2$Eout[2], y2)

        xs <- c(seg1_x, arc1_x[-1], seg2_x[-1], arc2_x[-1], seg3_x[-1])
        ys <- c(seg1_y, arc1_y[-1], seg2_y[-1], arc2_y[-1], seg3_y[-1])
      }

      df <- data.frame(x = xs, y = ys, PANEL = data$PANEL[i], group = grp)
      for (cl in style_cols) df[[cl]] <- style_vals[[cl]]
      df
    })

    do.call(rbind, out_list)
  }
)
