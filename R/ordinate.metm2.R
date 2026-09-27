#' Calculate ellipse coordinates for ordination plots
#'
#' Internal helper used by [ordinate.metm2()] to calculate confidence
#' or data ellipses for groups with at least three observations.
#'
#' @param dat A data frame containing `x`, `y`, and `Group` columns.
#' @param level Confidence or probability level.
#' @param type Ellipse type. One of `"confidence"` or `"data"`.
#' @param scale Numeric scaling factor applied to ellipse size.
#' @param npoint Number of points used to draw each ellipse.
#'
#' @return A data frame containing ellipse coordinates.
#'
#' @keywords internal
#' @noRd
.ord_ellipse <- function(
    dat,
    level = 0.95,
    type = c("confidence", "data"),
    scale = 1,
    npoint = 160
) {

  type <- match.arg(type)

  dplyr::bind_rows(
    lapply(
      split(dat, dat$Group),
      function(d) {

        z <- as.matrix(
          d[
            stats::complete.cases(
              d[, c("x", "y"), drop = FALSE]
            ),
            c("x", "y"),
            drop = FALSE
          ]
        )

        n <- nrow(z)

        ## 至少 3 个样本才计算椭圆
        if (n < 3) {
          return(NULL)
        }

        mu <- colMeans(z)

        cv <- stats::cov(z)

        ## 防止完全共线或零方差导致 eigen 出错
        if (
          any(!is.finite(cv)) ||
          all(abs(cv) < .Machine$double.eps)
        ) {
          return(NULL)
        }

        eg <- eigen(
          cv,
          symmetric = TRUE
        )

        eg$values <- pmax(
          eg$values,
          max(eg$values, 1) * 1e-12
        )

        k <- if (type == "confidence") {

          2 * (n - 1) /
            (n * (n - 2)) *
            stats::qf(
              level,
              2,
              n - 2
            )

        } else {

          stats::qchisq(
            level,
            df = 2
          )
        }

        A <- eg$vectors %*%
          diag(
            sqrt(
              eg$values * k
            ) * scale,
            2
          )

        th <- seq(
          0,
          2 * pi,
          length.out = npoint
        )

        xy <- t(
          matrix(
            mu,
            2,
            npoint
          ) +
            A %*%
            rbind(
              cos(th),
              sin(th)
            )
        )

        data.frame(
          x = xy[, 1],
          y = xy[, 2],
          Group = as.character(d$Group[1]),
          stringsAsFactors = FALSE
        )
      }
    )
  )
}


#' Ordination analysis with group ellipses
#'
#' Extends [ordinate.metm()] by adding confidence or data ellipses
#' to ordination plots. Ellipses are calculated only for groups
#' containing at least three valid observations.
#'
#' @param ... Arguments passed to [ordinate.metm()].
#' @param ellipse Logical. Whether to add group ellipses.
#' @param ellipse.level Numeric. Confidence or probability level used
#'   for ellipse construction. Default is `0.95`.
#' @param ellipse.type Character. Ellipse type. `"confidence"` draws
#'   a confidence ellipse around the group centroid, whereas `"data"`
#'   describes the dispersion of observations.
#' @param ellipse.scale Numeric. Multiplicative scaling factor for
#'   ellipse size.
#' @param ellipse.alpha Numeric. Transparency of the ellipse polygon.
#'
#' @return The same list returned by [ordinate.metm()], with ellipse
#'   layers added to available ggplot objects. An additional element
#'   named `ellipse` contains the calculated ellipse coordinates.
#'
#' @details
#' Groups with fewer than three valid samples are not assigned an ellipse.
#' The original ordination results from [ordinate.metm()] are otherwise
#' retained.
#'
#' @examples
#' \dontrun{
#' res <- ordinate.metm2(
#'   ps = ps,
#'   group = "Group",
#'   ellipse = TRUE,
#'   ellipse.level = 0.95,
#'   ellipse.type = "confidence"
#' )
#'
#' res[[1]]
#' }
#'
#' @export
ordinate.metm2 <- function(
    ...,
    ellipse = TRUE,
    ellipse.level = 0.95,
    ellipse.type = c("confidence", "data"),
    ellipse.scale = 1,
    ellipse.alpha = 0.12
) {

  ellipse.type <- match.arg(
    ellipse.type
  )

  ## 关键：直接调用同一 package 中的函数
  ## 不要写 EasyMultiOmics::ordinate.metm()
  res <- ordinate.metm(...)

  if (
    length(res) < 2 ||
    is.null(res[[2]])
  ) {
    return(res)
  }

  dat <- as.data.frame(
    res[[2]]
  )

  if (
    !ellipse ||
    !all(
      c(
        "x",
        "y",
        "Group"
      ) %in% names(dat)
    )
  ) {
    return(res)
  }

  ell <- .ord_ellipse(
    dat = dat,
    level = ellipse.level,
    type = ellipse.type,
    scale = ellipse.scale
  )

  if (
    is.null(ell) ||
    !nrow(ell)
  ) {
    res$ellipse <- data.frame()
    return(res)
  }

  add_ellipse <- function(p) {

    if (
      !inherits(
        p,
        "ggplot"
      )
    ) {
      return(p)
    }

    p +
      ggplot2::geom_polygon(
        data = ell,
        ggplot2::aes(
          x = x,
          y = y,
          group = Group,
          fill = Group
        ),
        inherit.aes = FALSE,
        alpha = ellipse.alpha,
        color = NA
      ) +
      ggplot2::geom_path(
        data = ell,
        ggplot2::aes(
          x = x,
          y = y,
          group = Group,
          color = Group
        ),
        inherit.aes = FALSE,
        linewidth = 0.7,
        linetype = 2
      )
  }

  ## 第一张图
  if (
    length(res) >= 1 &&
    inherits(
      res[[1]],
      "ggplot"
    )
  ) {
    res[[1]] <- add_ellipse(
      res[[1]]
    )
  }

  ## 如果第三个返回对象也是 ggplot，也添加椭圆
  if (
    length(res) >= 3 &&
    inherits(
      res[[3]],
      "ggplot"
    )
  ) {
    res[[3]] <- add_ellipse(
      res[[3]]
    )
  }

  res$ellipse <- ell

  res
}
