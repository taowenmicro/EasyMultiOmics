# feast_core.R
# 所有 FEAST 相关函数（对外接口 + 内部引擎 + 工具函数）

#' FEAST Microbial traceability analysis (package API)
#'
#' @description
#' 面向 R 包的 FEAST 包装函数。支持两种输入：
#' 1) 直接传 phyloseq 对象 `ps`（推荐）；或
#' 2) 传 `otu` / `map`（行=OTU/ASV，列=样本；map 行名=样本）。
#'
#' 函数内部会调用 FEAST 引擎（\code{FEAST_core}）对每个 sink 重复进行源追踪估计，
#' 并返回 source × sink 的比例表。
#'
#' @param ps A \code{phyloseq} object (可选).
#' @param otu OTU/ASV abundance table (taxa as rows, samples as columns).
#' @param map Sample metadata (rownames must be sample IDs).
#' @param group Column name in metadata used as group label.
#' @param sinkG The sink group label (single string).
#' @param sourceG The source group label(s), character vector.
#' @param rep Integer; biological replicates per group (auto-detected if NULL).
#' @param em_iter Integer; EM iterations for FEAST (default 1000).
#'
#' @return A data.frame: rows = sources(+unknown), cols = sink samples (replicates).
#' @examples
#' \dontrun{
#'   res <-  FEAST.micro(ps = ps, group = "Group",
#'                      sinkG = "WT", sourceG = c("KO","OE"))
#' }
#' @export
#' @import phyloseq vegan dplyr tibble magrittr
FEAST.micro <- function(ps = NULL,
                        otu = NULL,
                        map = NULL,
                        group = "Group",
                        sinkG = NULL,
                        sourceG = NULL,
                        rep = NULL,
                        em_iter = 1000) {

  if (is.null(ps)) {
    # --- 从 otu/map 组装 phyloseq ---
    stopifnot(!is.null(otu), !is.null(map))
    otu <- as.matrix(otu)
    if (is.null(rownames(otu)) || is.null(colnames(otu)))
      stop("otu must have rownames (taxa) and colnames (samples).")
    if (is.null(rownames(map)))
      stop("map must have rownames as sample IDs.")
    common <- intersect(colnames(otu), rownames(map))
    if (length(common) == 0)
      stop("No overlapping sample IDs between otu columns and map rownames.")
    otu <- otu[, common, drop = FALSE]
    map <- map[common, , drop = FALSE]
    ps <- phyloseq(otu_table(otu, taxa_are_rows = TRUE),
                   sample_data(as.data.frame(map)))
  } else {
    if (!inherits(ps, "phyloseq"))
      stop("ps must be a phyloseq object.")
  }

  # # --- 统一分组列 ---
  # md <- as.data.frame(sample_data(ps))
  # if (!group %in% colnames(md))
  #   stop(sprintf("group '%s' not found in sample_data(ps).", group))
  # md$Group <- as.factor(md[[group]])
  # sample_data(ps) <- md
  # --- 统一分组列 ---
  md <- phyloseq::sample_data(ps)
  # 转成纯 data.frame，去掉 <sample_data> 类
  md <- data.frame(md, check.names = FALSE, stringsAsFactors = FALSE)

  if (!group %in% colnames(md))
    stop(sprintf("group '%s' not found in sample_data(ps).", group))

  md$Group <- as.factor(md[[group]])

  # 把行名带入 id 列，然后用 base::order 排序，避免 dplyr::arrange 触发 vec_slice
  md$id <- rownames(md)
  md <- md[order(md$Group), , drop = FALSE]

  # 再写回 phyloseq
  sample_data(ps) <- md





  if (is.null(sinkG) || is.null(sourceG))
    stop("Please provide sinkG and sourceG.")
  sourceG <- as.character(sourceG)

  # --- OTU 矩阵 ---
  # 使用 vegan_otu 保持和原始代码一致
  otus <- as.data.frame(t(vegan_otu(ps)))
  otus <- t(as.matrix(otus))  # 基本等价于 vegan_otu(ps)

  # --- 样本顺序与分组 ---
  md$id <- rownames(md)
  md <- dplyr::arrange(md, .data$Group)
  envs <- md$id
  mu <- md$id[md$Group == sinkG]

  # --- 推断每组重复数 ---
  if (is.null(rep)) {
    rep <- length(md$Group) / length(unique(md$Group))
    rep <- as.integer(rep)
  }

  EM_iterations <- em_iter
  Proportions_est <- vector("list", length = rep)

  for (it in seq_len(rep)) {

    train.ix <- which(md$Group %in% sourceG &
                        md$id %in% md$id[seq(1, length(md$Group), rep) + (it - 1)])
    test.ix  <- which(md$Group == sinkG & md$id == mu[it])

    num_sources <- length(train.ix)
    COVERAGE <- min(rowSums(otus[c(train.ix, test.ix), , drop = FALSE]))

    sources <- as.matrix(vegan::rrarefy(otus[train.ix, , drop = FALSE], COVERAGE))
    sinks   <- as.matrix(vegan::rrarefy(t(as.matrix(otus[test.ix, , drop = FALSE])), COVERAGE))

    if (length(sourceG) == 1) {
      tmp <- as.data.frame(sources)
      sources <- rbind(tmp, tmp)
      rownames(sources) <- paste0(sourceG, 1:2)
      sources <- as.matrix(sources)
    }

    FEAST_out <- FEAST_core(source = sources,
                            sinks  = sinks,
                            env    = envs[train.ix],
                            em_itr = EM_iterations,
                            COVERAGE = COVERAGE)

    prop <- FEAST_out$data_prop[, 1]
    names(prop) <- c(as.character(envs[train.ix]), "unknown")

    if (length(prop) < num_sources + 1L) {
      tmp <- prop
      prop[num_sources]   <- NA_real_
      prop[num_sources+1] <- tmp[num_sources]
    }
    Proportions_est[[it]] <- prop
  }

  res <- as.data.frame(Proportions_est, check.names = FALSE)
  colnames(res) <- mu
  res
}

## ===========================
## 以下为 FEAST 内部引擎与工具函数
## 全部标记为 internal
## ===========================

#' @keywords internal
#' @noRd
change_C <- function(newcov, X){
  X <- t(as.matrix(X))
  idx <- 1:ncol(X)

  if (sum(X) > newcov) {
    while (sum(X) > newcov) {
      greaterone <- X > 1
      samps <- 20
      if (samps > length(X[greaterone]))
        samps <- length(X[greaterone])
      changeidx <- sample(idx[greaterone], samps, replace = FALSE)
      X[changeidx] <- X[changeidx] - 1
    }
  }

  if (sum(X) < newcov) {
    while (sum(X) < newcov) {
      greaterone <- X > 1
      samps <- 100
      if (samps > length(X[greaterone]))
        samps <- length(X[greaterone])
      changeidx <- sample(idx[greaterone], samps, replace = FALSE)
      X[changeidx] <- X[changeidx] + 1
    }
  }

  X
}

#' @keywords internal
#' @noRd
rarefy <- function(x, maxdepth){
  if (is.null(maxdepth)) return(x)

  if (!is.element(class(x)[1], c("matrix", "data.frame", "array")))
    x <- matrix(x, nrow = nrow(x))
  nr <- nrow(x)
  nc <- ncol(x)

  for(i in seq_len(nr)){
    if (sum(x[i, ]) > maxdepth) {
      prev.warn <- options()$warn
      options(warn = -1)
      s <- sample(nc, size = maxdepth, prob = x[i, ], replace = TRUE)
      options(warn = prev.warn)
      x[i, ] <- hist(s, breaks = seq(.5, nc + .5, 1), plot = FALSE)$counts
    }
  }
  x
}

#' @keywords internal
#' @noRd
jsdmatrix <- function(x){
  d <- matrix(0, nrow = nrow(x), ncol = nrow(x))
  for(i in 1:(nrow(x)-1)){
    for(j in (i+1):nrow(x)){
      d[i,j] <- jsd(x[i,], x[j,])
      d[j,i] <- d[i,j]
    }
  }
  d
}

#' @keywords internal
#' @noRd
jsd <- function(p, q){
  m <- (p + q) / 2
  (kld(p, m) + kld(q, m)) / 2
}

#' @keywords internal
#' @noRd
h <- function(x) { y <- x[x > 0]; -sum(y * log(y)) }

#' @keywords internal
#' @noRd
mult_JSD <- function(p, q) { h(q %*% p) - q %*% apply(p, 1, h) }

#' @keywords internal
#' @noRd
retrands <- function(V){
  unlist(lapply(c(V), function(x) runif(1, x + 1e-12, x + 1e-09)))
}

#' @keywords internal
#' @noRd
getR2 <- function(x, y){
  (cor(x, y))^2
}

#' @keywords internal
#' @noRd
E <- function(alphas, sources){
  nums <- (sapply(1:length(alphas),
                  function(n) Reduce("+", crossprod(as.numeric(alphas[n]), as.numeric(sources[[n]])))))
  denom <- (Reduce("+", nums))
  nums / denom
}

#' @keywords internal
#' @noRd
A <- function(alph, XO, raos){
  tmp <- crossprod(alph, XO/raos)
  tmp <- rapply(list(tmp), f = function(x) ifelse(is.nan(x), 0, x), how = "replace")
  tmp <- Reduce("+", unlist(tmp))
  tmp
}

#' @keywords internal
#' @noRd
M <- function(alphas, sources, sink, observed){
  newalphs <- c()
  rel_sink <- sink / sum(sink)

  if (sum(sources[[1]]) > 1) {
    sources <- lapply(sources, function(x) x / (sum(colSums(x))))
  }

  LOs <- lapply(sources, schur, b = rel_sink)
  BOs <- t(mapply(crossprod, x = sources, y = alphas))
  BOs <- split(BOs, seq(nrow(BOs)))
  BOs <- lapply(BOs, as.matrix)
  BOs <- lapply(BOs, t)
  num_list <- list()
  source_new <- list()

  for (i in 1:length(sources)) {
    num <- crossprod(alphas[i], (LOs[[i]] / (Reduce("+", BOs))))
    num <- rapply(list(num), f = function(x) ifelse(is.nan(x), 0, x), how = "replace")
    num_list[[i]] <- num[[1]][1, ] + observed[[i]][1, ]
    denom <- Reduce("+", unlist(num_list[[i]]))
    source_new[[i]] <- num_list[[i]] / denom
    source_new[[i]][is.na(source_new[[i]])] <- 0
  }

  sources <- source_new
  sources <- lapply(sources, t)
  XOs <- lapply(sources, schur, b = rel_sink)
  AOs <- t(mapply(crossprod, x = sources, y = alphas))
  AOs <- split(AOs, seq(nrow(AOs)))
  AOs <- lapply(AOs, as.matrix)
  AOs <- lapply(AOs, t)
  newAs <- c()

  for (i in 1:length(sources)) {
    newA <- crossprod(alphas[i], (XOs[[i]]/(Reduce("+", AOs))))
    newA <- rapply(list(newA), f = function(x) ifelse(is.nan(x), 0, x), how = "replace")
    newA <- Reduce("+", unlist(newA))
    newAs <- c(newAs, newA)
  }
  tot <- sum(newAs)
  Results <- list(new_alpha = newAs/tot, new_sources = sources)
  Results
}

#' @keywords internal
#' @noRd
do_EM <- function(alphas, sources, observed, sink, iterations){
  curalphas <- alphas
  newalphas <- alphas
  m_guesses <- c(alphas[1])

  for (itr in 1:iterations){
    curalphas <- E(newalphas, sources)
    tmp <- M(alphas = curalphas, sources = sources, sink = sink, observed = observed)
    newalphas <- tmp$new_alpha
    sources <- tmp$new_sources

    m_guesses <- c(m_guesses, newalphas[1])
    if (abs(m_guesses[length(m_guesses)] - m_guesses[length(m_guesses)-1]) <= 10^-6)
      break
  }
  toret <- c(newalphas)
  list(toret = toret, sources = sources)
}

#' @keywords internal
#' @noRd
M_basic <- function(alphas, sources, sink){
  XOs <- lapply(sources, schur, b = sink)
  AOs <- t(mapply(crossprod, x = sources, y = alphas))
  AOs <- split(AOs, seq(nrow(AOs)))
  AOs <- lapply(AOs, as.matrix)
  AOs <- lapply(AOs, t)
  newAs <- c()
  for (i in 1:length(sources)) {
    newA <- crossprod(alphas[i], (XOs[[i]]/(Reduce("+", AOs))))
    newA <- rapply(list(newA), f = function(x) ifelse(is.nan(x), 0, x), how = "replace")
    newA <- Reduce("+", unlist(newA))
    newAs <- c(newAs, newA)
  }
  tot <- sum(newAs)
  newAs / tot
}

#' @keywords internal
#' @noRd
do_EM_basic <- function(alphas, sources, sink, iterations){
  curalphas <- alphas
  newalphas <- alphas
  m_guesses <- c(alphas[1])

  for (itr in 1:iterations) {
    curalphas <- E(newalphas, sources)
    newalphas <- M_basic(curalphas, sources, sink)
    m_guesses <- c(m_guesses, newalphas[1])
    if (abs(m_guesses[length(m_guesses)] - m_guesses[length(m_guesses)-1]) <= 10^-6)
      break
  }
  c(newalphas)
}

#' @keywords internal
#' @noRd
source_process_nounknown <- function(train, envs, rarefaction_depth = 1000){
  train <- as.matrix(train)
  if (sum(as.integer(train) != as.numeric(train)) > 0) {
    stop("Data must be integral.")
  }
  envs <- factor(envs)
  train.envs <- sort(unique(levels(envs)))

  if (!is.null(rarefaction_depth) && rarefaction_depth > 0)
    train <- rarefy(train, rarefaction_depth)

  X <- t(sapply(split(data.frame(train), envs), colSums))
  rownames(X) <- train.envs
  X <- t(as.matrix(X))
  X
}

#' @keywords internal
#' @noRd
read_pseudo_data <- function(dataset){
  path_to_data <- "../data/"
  if (dataset == "DA") {
    df <- read.table(paste0(path_to_data,"DA_99_T_d10000_date_nan.txt"), fill = NA)
  } else if (dataset == "DB") {
    df <- read.table(paste0(path_to_data,"DB_99_T_d10000_date_nan.txt"), fill = NA)
  } else if (dataset == "F4") {
    df <- read.table(paste0(path_to_data,"F4_99_T_d10000_date_nan.txt"), fill = NA)
  } else {
    df <- read.table(paste0(path_to_data,"M3_99_T_d10000_date_nan.txt"), fill = NA)
  }
  df[complete.cases(df), ]
}

#' @keywords internal
#' @noRd
create_m <- function(num_sources, n, EPSILON){
  if (n == 1) {
    index <- sample(1:num_sources, 1)
    m_1 <- runif(min = 0.6, max = 0.9, n = 1)
    resid <- 1 - m_1
    other_ms <- resid / (num_sources - 1)
    m <- rep(NA_real_, num_sources)
    m[index] <- m_1
    m[is.na(m)] <- other_ms
  }

  if (n == 2) {
    index <- sample(1:num_sources, 2)
    m_1 <- runif(min = 0.1, max = 0.2, n = 1)
    m_2 <- runif(min = 0.4, max = 0.5, n = 1)
    resid <- 1 - (m_1 + m_2)
    other_ms <- resid / (num_sources - 2)
    m <- rep(NA_real_, num_sources)
    m[index] <- c(m_1, m_2)
    m[is.na(m)] <- other_ms
  }

  if (n == 3) {
    index <- sample(1:num_sources, 3)
    m_1 <- runif(min = 0.1, max = 0.5, n = 1)
    m_2 <- runif(min = 0.2, max = 0.25, n = 1)
    m_3 <- runif(min = 0.1, max = 0.15, n = 1)
    resid <- 1 - (m_1 + m_2 + m_3)
    other_ms <- runif(min = 0.001, max = resid/(num_sources-3), n = (num_sources-3))
    m <- rep(NA_real_, num_sources)
    m[index] <- c(m_1, m_2, m_3)
    m[is.na(m)] <- other_ms
    m <- m / sum(m)
  }

  subsum <- 0
  while ((subsum + 0.001) < EPSILON){
    tosub <- EPSILON - subsum
    tosub <- tosub / (num_sources + 1)
    mask <- m > tosub
    m[mask] <- m[mask] - tosub
    subsum <- subsum + length(m[mask]) * tosub
  }
  m <- c(m, EPSILON)
  m
}

#' @keywords internal
#' @noRd
unknown_initialize <- function(sources, sink, n_sources){
  unknown_source <- rep(0, length(sink))
  sum_sources <- apply(sources, 2, sum)
  for (j in seq_along(sum_sources)) {
    unknown_source[j] <- max(sink[j] - sum_sources[j], 0)
  }
  unknown_source
}

#' @keywords internal
#' @noRd
unknown_initialize_1 <- function(sources, sink, n_sources){
  unknown_source <- rep(0, length(sink))
  sources_sum <- apply(sources, 2, sum)

  ind_cor <- list()
  counter <- matrix(0, ncol = ncol(sources), nrow = nrow(sources))

  for (j in 1:n_sources){
    ind_cor[[j]] <- which(sources[j,] > 0)
    for (k in 1:ncol(sources)){
      if (sources[j,k] > 0){
        counter[j,k] <- counter[j,k] + 1
      }
    }
  }

  OTU_present_absent <- apply(counter, 2, sum)
  ind_cor_all <- which(OTU_present_absent >= round(n_sources * 0.8))

  if (length(ind_cor_all) > 1){
    cor_abundance <- round(apply(sources[, ind_cor_all], 2, median) / 2)
    unknown_source[ind_cor_all] <- cor_abundance
  }

  ind_no_known_source_abun <- which(sources_sum == 0)
  for (j in seq_along(ind_no_known_source_abun)){
    idx <- ind_no_known_source_abun[j]
    unknown_source[idx] <- max((sink[idx] - rpois(n = 1, lambda = 0.5)), 0)
  }
  unknown_source
}

#' @keywords internal
#' @noRd
unknown__initialize_1 <- function(sources, sink, n_sources){
  unknown_source <- rep(0, length(sink))
  sources_sum <- apply(sources, 2, sum)

  ind_cor <- list()
  counter <- matrix(0, ncol = ncol(sources), nrow = nrow(sources))

  for (j in 1:n_sources){
    ind_cor[[j]] <- which(sources[j,] > 0)
    for (k in 1:ncol(sources)){
      if (sources[j,k] > 0){
        counter[j,k] <- counter[j,k] + 1
      }
    }
  }

  OTU_present_absent <- apply(counter, 2, sum)
  ind_cor_all <- which(OTU_present_absent >= round(n_sources * 0.8))

  if (length(ind_cor_all) > 1){
    cor_abundance <- apply(sources[, ind_cor_all], 2, median)
    unknown_source[ind_cor_all] <- cor_abundance
  }

  ind_no_known_source_abun <- which(sources_sum == 0)
  for (j in seq_along(ind_no_known_source_abun)){
    idx <- ind_no_known_source_abun[j]
    unknown_source[idx] <- max(round(sink[idx] + rnorm(n = 1)), 0)
  }
  unknown_source
}

#' Internal FEAST engine
#' @keywords internal
#' @noRd
FEAST_core <- function(source,
                       sinks,
                       em_itr = 1000,
                       env = rownames(source),
                       include_epsilon = TRUE,
                       COVERAGE,
                       unknown_initialize = 0){

  tmp <- source
  test_zeros <- apply(tmp, 1, sum)
  ind_to_use <- which(test_zeros > 0)
  ind_zero   <- which(test_zeros == 0)

  source <- tmp[ind_to_use, , drop = FALSE]

  totalsource <- source
  totalsource <- as.matrix(totalsource)
  sources <- split(totalsource, seq(nrow(totalsource)))
  sources <- lapply(sources, as.matrix)
  dists <- lapply(sources, function(x) x/(sum(colSums(x))))
  totaldist <- t(Reduce("cbind", dists))
  sinks <- matrix(sinks, nrow = 1, ncol = ncol(totalsource))

  num_sources <- nrow(source)
  envs_simulation <- seq_len(num_sources)

  source_old <- source
  totalsource_old <- totalsource

  source_old <- lapply(source_old, t)
  source_old <- split(totalsource_old, seq(nrow(totalsource_old)))
  source_old <- lapply(source_old, as.matrix)

  if (include_epsilon) {

    source_2 <- list()
    totalsource_2 <- matrix(NA_real_,
                            ncol = ncol(totalsource_old),
                            nrow = nrow(totalsource_old) + 1L)

    for (j in 1:num_sources){
      source_2[[j]] <- source_old[[j]]
      totalsource_2[j, ] <- totalsource_old[j, ]
    }

    sinks_rarefy <- rarefy(matrix(sinks, nrow = 1),
                           maxdepth = apply(totalsource_old, 1, sum)[1])

    if (unknown_initialize == 1)
      unknown_source_1 <- unknown_initialize_1(
        sources = totalsource[1:num_sources, , drop = FALSE],
        sink    = as.numeric(sinks),
        n_sources = num_sources
      )

    if (unknown_initialize == 0)
      unknown_source_1 <- unknown_initialize(
        sources = totalsource[1:num_sources, , drop = FALSE],
        sink    = as.numeric(sinks),
        n_sources = num_sources
      )

    unknown_source <- unknown_source_1 + rpois(n = length(sinks), lambda = 0.5)

    unknown_source_rarefy <- rarefy(matrix(unknown_source, nrow = 1),
                                    maxdepth = COVERAGE)
    source_2[[num_sources + 1L]] <- t(unknown_source_rarefy)
    totalsource_2[num_sources + 1L, ] <- t(unknown_source_rarefy)
    totalsource <- totalsource_2

    source <- lapply(source_2, t)
    source <- split(totalsource, seq(nrow(totalsource_2)))
    source <- lapply(source_2, as.matrix)

    envs_simulation <- seq_len(num_sources + 1L)
  }

  samps <- source
  samps <- lapply(samps, t)
  observed_samps <- samps
  observed_samps[[(num_sources + 1L)]] <- t(rep(0, ncol(samps[[1]])))

  initalphs <- runif(num_sources + 1L, 0.0, 1.0)
  initalphs <- initalphs / sum(initalphs)

  sink_em <- as.matrix(sinks)
  pred_em <- do_EM_basic(alphas = initalphs,
                         sources = samps,
                         sink = sink_em,
                         iterations = em_itr)

  tmp2 <- do_EM(alphas = initalphs,
                sources = samps,
                sink = sink_em,
                iterations = em_itr,
                observed = observed_samps)
  pred_emnoise <- tmp2$toret

  k <- 1L
  pred_emnoise_all <- c()
  pred_em_all <- c()

  for (j in seq_along(env)){
    if (j %in% ind_to_use) {
      pred_emnoise_all[j] <- pred_emnoise[k]
      pred_em_all[j]      <- pred_em[k]
      k <- k + 1L
    } else {
      pred_emnoise_all[j] <- 0
      pred_em_all[j]      <- 0
    }
  }

  pred_emnoise_all[length(env) + 1L] <- pred_emnoise[k]
  pred_em_all[length(env) + 1L]      <- pred_em[k]

  names(pred_emnoise_all) <- c(env, "unknown")
  names(pred_em_all)      <- c(env, "unknown")

  Results <- list(
    unknown_source       = unknown_source,
    unknown_source_rarefy = unknown_source_rarefy,
    data_prop            = data.frame(pred_emnoise_all, pred_em_all)
  )
  Results
}
