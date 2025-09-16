#' @export
diff_all_methods <- function(ps, group = "Group", alpha = 0.05) {
  results <- list()

  # ---- 1. ALDEx2 ----
  results[["ALDEx2"]] <- tryCatch({
    aldex2.micro(ps, group, alpha)$diff.tab %>%
      dplyr::select(micro, adjust.p, method)
  }, error=function(e) NULL)

  # ---- 2. ANCOM fast ----
  results[["ANCOMII.fast"]] <- tryCatch({
    ancom.micro.fast(ps, group, alpha)$diff.tab %>%
      dplyr::select(micro, adjust.p, method)
  }, error=function(e) NULL)

  # ---- 3. corncob ----
  results[["corncob"]] <- tryCatch({
    df <- corncob.micro(ps, group, alpha)$diff.tab
    colnames(df) <- c("micro", "adjust.p", "method")
    df
  }, error=function(e) NULL)

  # ---- 4. edgeR ----
  results[["edgeR"]] <- tryCatch({
    edgeR.micro(ps, group, alpha)$diff.tab %>%
      dplyr::select(micro, adjust.p, method)
  }, error=function(e) NULL)

  # ---- 5. DESeq2 ----
  results[["DESeq2"]] <- tryCatch({
    deseq2.micro(ps, group, alpha)$diff.tab %>%
      dplyr::select(micro, adjust.p, method)
  }, error=function(e) NULL)

  # ---- 6. LEfSe ----
  results[["LEfSe"]] <- tryCatch({
    lefse.micro(ps, group, alpha)$diff.tab %>%
      dplyr::select(micro, adjust.p, method)
  }, error=function(e) NULL)

  # ---- 7. limma voom (TMM) ----
  results[["limma.voom.TMM"]] <- tryCatch({
    limma.v.TMM.micro(ps, group, alpha, method="TMM")$diff.tab %>%
      dplyr::select(micro, adjust.p, method)
  }, error=function(e) NULL)

  # ---- 8. limma voom (TMMwsp) ----
  results[["limma.voom.TMMwsp"]] <- tryCatch({
    limma.v.TMM.micro(ps, group, alpha, method="TMMwsp")$diff.tab %>%
      dplyr::select(micro, adjust.p, method)
  }, error=function(e) NULL)

  # ---- 9. Maaslin2 ----
  results[["Maaslin2"]] <- tryCatch({
    maaslin2.micro(ps, group, alpha, rare=FALSE)$diff.tab %>%
      dplyr::select(micro, adjust.p, method)
  }, error=function(e) NULL)

  # ---- 10. Maaslin2 (稀释) ----
  results[["Maaslin2.rare"]] <- tryCatch({
    maaslin2.micro(ps, group, alpha, rare=TRUE)$diff.tab %>%
      dplyr::select(micro, adjust.p, method)
  }, error=function(e) NULL)

  # ---- 11. metagenomeSeq ----
  results[["metagenomeSeq"]] <- tryCatch({
    metaSeq.micro(ps, group, alpha)$diff.tab %>%
      dplyr::select(micro, adjust.p, method)
  }, error=function(e) NULL)

  # ---- 12. t.test (稀释) ----
  results[["t.test.rare"]] <- tryCatch({
    t.test.micro(ps, group, alpha)$diff.tab %>%
      dplyr::select(micro, adjust.p, method)
  }, error=function(e) NULL)

  # ---- 13. wilcox.test (稀释) ----
  results[["wilcox.rare"]] <- tryCatch({
    wilcox.sampl.micro(ps, group, alpha)$diff.tab %>%
      dplyr::select(micro, adjust.p, method)
  }, error=function(e) NULL)

  # ---- 14. wilcox.test (CLR) ----
  results[["wilcox.CLR"]] <- tryCatch({
    wilcox.clr.micro(ps, group, alpha)$diff.tab %>%
      dplyr::select(micro, adjust.p, method)
  }, error=function(e) NULL)

  # ---- 15. ANCOMBC2 ----
  results[["ANCOMBC2"]] <- tryCatch({
    df <- ancombc2.micro(ps, group, alpha)$diff.tab
    df$contrast <- NULL
    df %>% dplyr::select(micro, adjust.p, method)
  }, error=function(e) NULL)

  # ---- 16. ZicoSeq ----
  results[["ZicoSeq"]] <- tryCatch({
    df <- zicoseq.micro(ps, group, alpha)$diff.tab
    df$p.fwer <- NULL; df$p.raw <- NULL
    df %>% dplyr::select(micro, adjust.p, method)
  }, error=function(e) NULL)

  # ---- 17. SongbirdR ----
  results[["songbirdR"]] <- tryCatch({
    df <- songbirdR.micro(ps, group, collapse=TRUE, pseudo_p=TRUE, sig_alpha=alpha)$diff.tab
    df <- df %>% dplyr::filter(significant == TRUE)
    df$Group <- NULL; df$coef <- NULL; df$significant <- NULL
    df %>% dplyr::select(micro, adjust.p, method)
  }, error=function(e) NULL)

  #
  all_tabs <- dplyr::bind_rows(results, .id = "method_source") %>%
    dplyr::filter(!is.na(micro) & micro != "")

  # 统计每个 OTU 被多少方法识别
  summary_tab <- all_tabs %>%
    dplyr::group_by(micro) %>%
    dplyr::summarise(
      n_methods = dplyr::n_distinct(method),
      methods = paste(unique(method), collapse = "; "),
      min_p = suppressWarnings(min(adjust.p, na.rm = TRUE))
    ) %>%
    dplyr::arrange(desc(n_methods), min_p)

  return(list(
    summary = summary_tab,
    details = all_tabs,
    by_method = results
  ))
}







# --- 1. ALDEx2 ---
#' @export
aldex2.micro <- function(ps, group = "Group", alpha = 0.05, test = "t") {
  ASV_table <- ps %>% filter_taxa(function(x) sum(x) > 0 , TRUE) %>%
    vegan_otu() %>% t() %>% as.data.frame()
  groupings <- sample_data(ps)
  results <- ALDEx2::aldex(reads = ASV_table,
                           conditions = groupings[[group]] %>% as.vector(),
                           mc.samples = 128, test = test,
                           effect = TRUE, denom = "all")
  tab.d <- results %>% as.data.frame() %>%
    dplyr::filter(we.ep < alpha) %>% tibble::rownames_to_column("id") %>%
    dplyr::select(id,we.ep) %>% dplyr::rename(OTU = id, p = we.ep) %>%
    dplyr::mutate(group = "ALDEx2")
  list(tab.Aldex2 = results,
       diff.tab = data.frame(micro = tab.d$OTU, method = tab.d$group, adjust.p = tab.d$p))
}

# --- 2. ANCOM-II ---
#' @export
ancom.micro <- function(ps, group = "Group", alpha = 0.05,
                        p_adj_method = "BH", ANCOM.value = "detected_0.6") {
  ASV_table <- vegan_otu(ps) %>% t() %>% as.data.frame()
  groupings <- sample_data(ps); groupings$Sample <- rownames(groupings)
  prepro <- feature_table_pre_process(feature_table = ASV_table,
                                      meta_data = groupings, sample_var = 'Sample',
                                      group_var = NULL, out_cut = 0.05, zero_cut = 0.90,
                                      lib_cut = 1000, neg_lb=FALSE)
  res <- ANCOM(feature_table = prepro$feature_table,
               meta_data = prepro$meta_data,
               struc_zero = prepro$structure_zeros,
               main_var = group, p_adj_method = p_adj_method, alpha = alpha)
  dat <- res$out
  tab.d <- dat %>% dplyr::filter(.data[[ANCOM.value]] == TRUE) %>%
    dplyr::select(taxa_id,.data[[ANCOM.value]]) %>%
    dplyr::rename(OTU = taxa_id, p = !!ANCOM.value) %>% dplyr::mutate(group = "ANCOMII")
  list(tab.ANCOM = dat,
       diff.tab = data.frame(micro = tab.d$OTU, method = tab.d$group, adjust.p = tab.d$p))
}

# --- 3. corncob ---
#' @export
corncob.micro <- function(ps, group = "Group", alpha = 0.05) {
  my_formula <- as.formula(paste("~", group))

  results <- corncob::differentialTest(
    formula = my_formula,
    phi.formula = my_formula,
    phi.formula_null = my_formula,
    formula_null = ~1,
    test = "Wald",
    data = ps,
    boot = FALSE,
    fdr_cutoff = alpha
  )

  dat <- data.frame(
    OTU = names(results$p_fdr),
    p_fdr = results$p_fdr
  )

  # 只保留显著差异的 OTU
  tab.d <- dat %>%
    dplyr::filter(p_fdr < alpha) %>%
    dplyr::mutate(
      p = p_fdr,
      group = "corncob"
    ) %>%
    dplyr::select(OTU, p, group)

  return(list(
    tab.corncob = dat,
    diff.tab = tab.d
  ))
}

# --- 4. edgeR ---

#' @export
edgeR.micro <- function(ps, group = "Group", alpha = 0.05) {

  phyloseq_to_edgeR <- function(physeq, group, method = "RLE") {
    if (!taxa_are_rows(physeq)) { physeq <- t(physeq) }
    x <- as(otu_table(physeq), "matrix") + 1

    # 处理 group
    if (length(group) == 1 && nsamples(physeq) > 1) {
      group <- get_variable(physeq, group)
    }

    # taxonomy 可选
    taxonomy <- NULL
    if (!is.null(phyloseq::tax_table(physeq))) {
      taxonomy <- as(phyloseq::tax_table(physeq), "matrix") %>% as.data.frame()
    }

    # DGEList
    y <- edgeR::DGEList(counts = x, group = group, genes = taxonomy)
    z <- edgeR::calcNormFactors(y, method = method)
    edgeR::estimateTagwiseDisp(edgeR::estimateCommonDisp(z))
  }

  test <- phyloseq_to_edgeR(ps, group)
  et <- edgeR::exactTest(test)
  tt <- edgeR::topTags(et, n = nrow(test$table), adjust.method = "fdr")
  res <- tt@.Data[[1]]

  tab.d <- res %>%
    tibble::rownames_to_column("OTU") %>%
    dplyr::filter(FDR < alpha) %>%
    dplyr::transmute(micro = OTU, method = "edgeR", adjust.p = FDR)

  return(list(
    tab.edgeR = res,
    diff.tab = tab.d
  ))
}

# --- 5. DESeq2 ---
#' @export
deseq2.micro <- function(ps, group = "Group", alpha = 0.05) {
  my_formula <- as.formula(paste("~",group))
  ASV_table <- vegan_otu(ps) %>% t() %>% as.data.frame()
  groupings <- sample_data(ps)
  dds <- DESeq2::DESeqDataSetFromMatrix(ASV_table, colData=groupings, design=my_formula)
  dds_res <- DESeq2::DESeq(dds, sfType="poscounts")
  res <- DESeq2::results(dds_res, tidy=TRUE, format="DataFrame")
  rownames(res) <- res$row; res <- res[,-1]
  tab.d <- res %>% tibble::rownames_to_column("id") %>%
    dplyr::filter(padj < alpha) %>% dplyr::rename(OTU = id, p = padj) %>%
    dplyr::mutate(group = "DESeq2")
  list(tab.DESeq2 = res,
       diff.tab = data.frame(micro = tab.d$OTU, method = tab.d$group, adjust.p = tab.d$p))
}

# --- 6. LEfSe ---
#' @export
lefse.micro <- function(ps, group = "Group", alpha = 0.05) {
  tablda <- LDA_Micro(ps, group = group, p.lvl = alpha, lda.lvl = 0,
                      seed = 11, adjust.p = FALSE)
  dat <- tablda[[2]]; dat$ID <- row.names(dat)
  dat2 <- dat %>% dplyr::filter(stringr::str_detect(ID,"st__"))
  dat2$OTU <- gsub("st__","",dat2$ID)
  print(1)
  tab.d <- dat2 %>% dplyr::filter(Pvalues < alpha) %>%
    dplyr::select(OTU,Pvalues) %>%
    dplyr::rename(p = Pvalues) %>%
    dplyr::mutate(group = "LEfSe")
  list(tab.LEfSe = dat2,
       diff.tab = data.frame(micro = tab.d$OTU, method = tab.d$group, adjust.p = tab.d$p))
}

# --- 7. limma voom ---
#' @export
limma.v.TMM.micro <- function(ps, group="Group", alpha=0.05, method="TMM") {
  # 构建公式
  my_formula <- as.formula(paste("~", paste(group, collapse=" + ")))

  # OTU 表
  ASV_table <- ps %>%
    phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) %>%
    vegan_otu() %>% t() %>% as.data.frame()

  # 提取 metadata，转成干净的 data.frame
  groupings <- as(sample_data(ps), "data.frame")

  # 确保因子型
  for (g in group) {
    if (!is.factor(groupings[[g]])) {
      groupings[[g]] <- factor(groupings[[g]])
    }
  }

  # 构建 DGEList
  DGE_LIST <- edgeR::DGEList(ASV_table)
  DGE_LIST_Norm <- edgeR::calcNormFactors(DGE_LIST, method=method)

  # 设计矩阵
  mm <- model.matrix(my_formula, data=groupings)

  # voom + limma
  voomvoom <- limma::voom(DGE_LIST_Norm, mm, plot=FALSE)
  fit <- limma::lmFit(voomvoom, mm)
  fit <- limma::eBayes(fit)

  res <- limma::topTable(fit, coef=2, n=nrow(DGE_LIST_Norm), sort.by="none")

  tab.d <- res %>%
    tibble::rownames_to_column("OTU") %>%
    dplyr::filter(adj.P.Val < alpha) %>%
    dplyr::transmute(micro = OTU,
                     method = paste0("limma.voom.", method),
                     adjust.p = adj.P.Val)

  return(list(
    tab.limma = res,
    diff.tab = tab.d
  ))
}

# --- 8. Maaslin2 ---
#' @export
maaslin2.micro <- function(ps, group = "Group", alpha = 0.05, rare = FALSE) {

  # OTU 表
  if (rare) {
    ASV_table <- ps %>%
      scale_micro(method = "sampling") %>%
      vegan_otu() %>% t()
    tem <- "Maaslin2.rare"
  } else {
    ASV_table <- ps %>%
      filter_taxa(function(x) sum(x) > 0, TRUE) %>%
      vegan_otu()  %>% t()
    tem <- "Maaslin2"
  }

  # Maaslin2 要求：行 = feature (OTU)，列 = sample
  ASV_table <- as.data.frame(ASV_table)

  # metadata
  groupings <- as(sample_data(ps), "data.frame")

  # 确保分组是因子
  if (!is.factor(groupings[[group]])) {
    groupings[[group]] <- factor(groupings[[group]])
  }

  # 保证样本对齐
  common.samples <- intersect(rownames(groupings), colnames(ASV_table))
  groupings <- groupings[common.samples, , drop = FALSE]
  ASV_table <- ASV_table[, common.samples, drop = FALSE]

  # Maaslin2 要求 input_data 行名必须是 OTU
  stopifnot(!is.null(rownames(ASV_table)))

  # 运行 Maaslin2
  fit_data <- Maaslin2::Maaslin2(
    input_data = ASV_table,
    input_metadata = groupings,
    output = tem,
    transform = "AST",
    fixed_effects = group,
    standardize = FALSE,
    plot_heatmap = FALSE,
    plot_scatter = FALSE
  )

  dat <- fit_data$results

  tab.d <- dat %>%
    dplyr::filter(qval < alpha) %>%
    dplyr::select(feature, qval) %>%
    dplyr::rename(OTU = feature, p = qval) %>%
    dplyr::mutate(group = tem)

  unlink(tem, recursive = TRUE)

  return(list(
    tab.Maaslin2 = dat,
    diff.tab = data.frame(micro = tab.d$OTU,
                          method = tab.d$group,
                          adjust.p = tab.d$p)
  ))
}

# --- 9. metagenomeSeq ---
#' @export
metaSeq.micro <- function(ps, group = "Group", alpha = 0.05) {
  # 提取 OTU 表，行 = OTU，列 = sample
  ASV_table <- ps %>%
    filter_taxa(function(x) sum(x) > 0, TRUE) %>%
    vegan_otu() %>% t()
  ASV_table <- as.data.frame(ASV_table)

  # metadata
  groupings <- as(sample_data(ps), "data.frame")

  # 保证样本名对齐
  common.samples <- intersect(rownames(groupings), colnames(ASV_table))
  groupings <- groupings[common.samples, , drop = FALSE]
  ASV_table <- ASV_table[, common.samples, drop = FALSE]

  # featureData 必须和 OTU 表行名一致
  feature_data <- AnnotatedDataFrame(data.frame(OTU = rownames(ASV_table)))
  rownames(feature_data) <- rownames(ASV_table)

  # phenoData 必须和 OTU 表列名一致
  pheno <- AnnotatedDataFrame(groupings)
  rownames(pheno) <- rownames(groupings)

  # 构建 MRexperiment
  test_obj <- metagenomeSeq::newMRexperiment(counts = ASV_table,
                                             phenoData = pheno,
                                             featureData = feature_data)

  # 归一化
  p <- metagenomeSeq::cumNormStat(test_obj)
  test_obj_norm <- metagenomeSeq::cumNorm(test_obj, p = p)

  # 建模
  mod <- model.matrix(as.formula(paste("~", group)), data = groupings)
  regres <- metagenomeSeq::fitFeatureModel(test_obj_norm, mod)

  res_table <- metagenomeSeq::MRfulltable(regres, number = nrow(ASV_table))

  tab.d <- res_table %>%
    tibble::rownames_to_column("OTU") %>%
    dplyr::filter(adjPvalues < alpha) %>%
    dplyr::transmute(micro = OTU, method = "metagenomeSeq", adjust.p = adjPvalues)

  return(list(
    tab.metaSeq = res_table,
    diff.tab = tab.d
  ))
}

# --- 10. t.test (rare) ---
#' @export
t.test.micro <- function(ps, group="Group", alpha=0.05) {
  ASV_table <- ps %>% scale_micro(method="sampling") %>% vegan_otu() %>% t() %>% as.data.frame()
  g <- sample_data(ps)[[group]]
  pvals <- apply(ASV_table,1,function(x) t.test(x~g,exact=FALSE)$p.value)
  dat <- data.frame(id=rownames(ASV_table),p=pvals)
  tab.d <- dat %>% dplyr::filter(p<alpha) %>% dplyr::rename(OTU=id) %>% dplyr::mutate(group="t.test.rare")
  list(tab.ttest=dat,
       diff.tab=data.frame(micro=tab.d$OTU,method=tab.d$group,adjust.p=tab.d$p))
}

# --- 11. Wilcoxon (rare) ---
#' @export
wilcox.sampl.micro <- function(ps, group="Group", alpha=0.05) {
  ASV_table <- ps %>% scale_micro(method="sampling") %>% vegan_otu() %>% t() %>% as.data.frame()
  g <- sample_data(ps)[[group]]
  pvals <- apply(ASV_table,1,function(x) wilcox.test(x~g,exact=FALSE)$p.value)
  dat <- data.frame(id=rownames(ASV_table),p=pvals)
  tab.d <- dat %>% dplyr::filter(p<alpha) %>% dplyr::rename(OTU=id) %>% dplyr::mutate(group="wilcox.rare")
  list(tab.wilcox=dat,
       diff.tab=data.frame(micro=tab.d$OTU,method=tab.d$group,adjust.p=tab.d$p))
}

# --- 12. Wilcoxon (CLR) ---
#' @export
wilcox.clr.micro <- function(ps, group="Group", alpha=0.05) {
  ASV_table <- vegan_otu(ps)%>%t()%>%as.data.frame()
  CLR_table <- data.frame(apply(ASV_table+1,2,function(x){log(x)-mean(log(x))}))
  g <- sample_data(ps)[[group]]
  pvals <- apply(CLR_table,1,function(x) wilcox.test(x~g,exact=FALSE)$p.value)
  dat <- data.frame(id=rownames(CLR_table),p=pvals)
  tab.d <- dat %>% dplyr::filter(p<alpha) %>% dplyr::rename(OTU=id) %>% dplyr::mutate(group="wilcox.CLR")
  list(tab.wilcoxCLR=dat,
       diff.tab=data.frame(micro=tab.d$OTU,method=tab.d$group,adjust.p=tab.d$p))
}

# --- 13. ANCOMBC ---
#' @export
ancombc2.micro <- function(ps, group = "Group", alpha = 0.05) {
  # 转换 phyloseq 为 TreeSummarizedExperiment
  message("正在将 phyloseq 转换为 TreeSummarizedExperiment...")
  tse <- mia::convertFromPhyloseq(ps)

  # 运行 ANCOM-BC2 分析
  message("正在运行 ANCOM-BC2 分析...")
  out <- ANCOMBC::ancombc2(
    data = tse,
    assay_name = "counts",
    tax_level = NULL,
    fix_formula = group,
    p_adj_method = "holm",
    prv_cut = 0.10,
    lib_cut = 1000,
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = alpha,
    global = TRUE,
    group = group,
    n_cl = 8
  )

  # 检查 out 是否包含 res 对象
  if (!"res" %in% names(out)) {
    stop("ANCOM-BC2 输出中未找到 'res' 对象。请检查 ANCOMBC::ancombc2 的输出。")
  }

  # 输出 res 的列名以便调试
  message("out$res 的列名: ", paste(colnames(out$res), collapse = ", "))

  # 检查 taxon 列是否存在
  if (!"taxon" %in% colnames(out$res)) {
    stop("out$res 中未找到 'taxon' 列。可用列名: ",
         paste(colnames(out$res), collapse = ", "))
  }

  # 检查是否存在以 q_ 开头的列
  q_cols <- colnames(out$res)[grepl("^q_", colnames(out$res))]
  if (length(q_cols) == 0) {
    warning("未找到以 'q_' 开头的列，可能影响结果过滤。")
  }

  # 处理差异丰度结果
  tab.d <- out$res %>%
    dplyr::select(taxon, dplyr::starts_with("q_")) %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with("q_"),
      names_to = "contrast",
      names_prefix = "q_",
      values_to = "q_val"
    ) %>%
    dplyr::filter(!is.na(q_val) & q_val < alpha) %>%
    dplyr::rename(OTU = taxon, p = q_val) %>%
    dplyr::mutate(method = "ancombc2")

  # 返回结果
  return(list(
    tab.ancombc2 = out,
    diff.tab = data.frame(
      micro = tab.d$OTU,
      contrast = tab.d$contrast,
      method = tab.d$method,
      adjust.p = tab.d$p
    )
  ))
}
# --- 14. ZicoSeq ---
#' @export
zicoseq.micro <- function(ps, group = "Group", alpha = 0.05) {
  comm <- ps %>%
    filter_taxa(function(x) sum(x) > 0, TRUE) %>%
    vegan_otu() %>% t()

  # 去掉全零和无差异的 feature
  comm <- comm[rowSums(comm) > 0, , drop = FALSE]
  comm <- comm[apply(comm, 1, var) > 0, , drop = FALSE]

  meta.dat <- as(sample_data(ps), "data.frame")

  ZicoSeq.obj <- GUniFrac::ZicoSeq(
    meta.dat = meta.dat,
    feature.dat = comm,
    grp.name = group,
    feature.dat.type = "count",
    is.winsor = TRUE,
    outlier.pct = 0.03,
    is.post.sample = TRUE,
    post.sample.no = 25,
    link.func = list(function(x) x^0.5),
    perm.no = 99,
    ref.pct = 0.5,
    stage.no = 6,
    excl.pct = 0.2,
    is.fwer = TRUE
  )

  dat <- data.frame(
    id = names(ZicoSeq.obj$p.raw),
    p.raw = ZicoSeq.obj$p.raw,
    p.adj = ZicoSeq.obj$p.adj.fdr,
    p.fwer = ZicoSeq.obj$p.adj.fwer
  )

  tab.d <- dat %>%
    dplyr::filter(p.adj < alpha) %>%
    dplyr::mutate(method = "ZicoSeq") %>%
    dplyr::rename(micro = id, adjust.p = p.adj)

  return(list(
    tab.ZicoSeq = dat,
    diff.tab = tab.d
  ))
}

#' @export
LDA_Micro = function(ps = ps,
                     group = "Group",
                     Top = 100,
                     p.lvl = 0.05,
                     lda.lvl = 2,
                     seed = 11,
                     adjust.p = F
){

  ps = ps %>%  filter_taxa(function(x) sum(x ) > 0 , TRUE)
  alltax = ps %>%
    ggClusterNet::filter_OTU_ps(Top) %>%
    ggClusterNet::vegan_tax() %>%
    as.data.frame()
  alltax$OTU = row.names(alltax)

  alltax$Kingdom = paste(alltax$Kingdom,sep = "_Rank_")
  alltax$Phylum = paste(alltax$Kingdom,alltax$Phylum,sep = "_Rank_")
  alltax$Class = paste(alltax$Phylum,alltax$Class,sep = "_Rank_")
  alltax$Order = paste(alltax$Class,alltax$Order,sep = "_Rank_")
  alltax$Family = paste(alltax$Order,alltax$Family,sep = "_Rank_")
  alltax$Genus = paste(alltax$Family,alltax$Genus,sep = "_Rank_")
  alltax$Species = paste(alltax$Genus,alltax$Species,sep = "_Rank_")


  otu = ps %>%
    ggClusterNet::filter_OTU_ps(Top) %>%
    ggClusterNet::vegan_otu() %>%
    t() %>%
    as.data.frame()

  otu_tax = merge(otu,alltax,by = "row.names",all = F)
  head(otu_tax)

  rank1 <- otu_tax %>%
    dplyr::group_by(Kingdom) %>%
    dplyr::summarise_if(is.numeric, sum, na.rm = TRUE)
  colnames(rank1)[1] = "id"
  rank1$id = paste("k__",rank1$id,sep = "")
  rank2 <- otu_tax %>%
    dplyr::group_by(Phylum) %>%
    dplyr::summarise_if(is.numeric, sum, na.rm = TRUE)
  colnames(rank2)[1] = "id"
  rank2$id = paste("p__",rank2$id,sep = "")
  rank3 <- otu_tax %>%
    dplyr::group_by(Class) %>%
    dplyr::summarise_if(is.numeric, sum, na.rm = TRUE)
  colnames(rank3)[1] = "id"
  rank3$id = paste("c__",rank3$id,sep = "")

  rank4 <- otu_tax %>%
    dplyr::group_by(Order) %>%
    dplyr::summarise_if(is.numeric, sum, na.rm = TRUE)
  colnames(rank4)[1] = "id"
  rank4$id = paste("o__",rank4$id,sep = "")

  rank5 <- otu_tax %>%
    dplyr::group_by(Family) %>%
    dplyr::summarise_if(is.numeric, sum, na.rm = TRUE)
  colnames(rank5)[1] = "id"
  rank5$id = paste("f__",rank5$id,sep = "")

  rank6 <- otu_tax %>%
    dplyr::group_by(Genus) %>%
    dplyr::summarise_if(is.numeric, sum, na.rm = TRUE)
  colnames(rank6)[1] = "id"
  rank6$id = paste("g__",rank6$id,sep = "")

  rank7 <- otu_tax %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise_if(is.numeric, sum, na.rm = TRUE)
  colnames(rank7)[1] = "id"
  rank7$id = paste("s__",rank7$id,sep = "")

  rank8 <- otu_tax %>%
    dplyr::group_by(OTU) %>%
    dplyr::summarise_if(is.numeric, sum, na.rm = TRUE)
  colnames(rank8)[1] = "id"
  rank8$id = paste("st__",rank8$id,sep = "")

  # 合并8个分类级
  all = rbind(rank1,rank2,rank3,rank4,rank5,rank6,rank7,rank8)
  head(all)


  #
  data1 = as.data.frame(all)
  row.names(data1) = data1$id
  data1$id = NULL

  #-构建phylose对象

  ps_G_graphlan = phyloseq::phyloseq(phyloseq::otu_table(as.matrix(data1),taxa_are_rows = TRUE),
                                     phyloseq::sample_data(ps)) %>%  filter_taxa(function(x) sum(x ) > 0, TRUE)
  ps_G_graphlan

  #----提取OTU表格

  otu = as.data.frame((ggClusterNet::vegan_otu(ps_G_graphlan)))
  otu[otu==0] <- 1
  otu = otu[ colMeans(otu) != 1]


  map = as.data.frame(phyloseq::sample_data(ps_G_graphlan))
  # otu = (otu_table)
  claslbl= map[,group] %>% as.vector() %>% .[[1]] %>% as.factor()
  # claslbl= map$Group %>% as.factor()
  set.seed(seed)
  #KW rank sum test

  rawpvalues <- apply(otu, 2, function(x) kruskal.test(x, claslbl)$p.value);
  #--得到计算后得到的p值
  ord.inx <- order(rawpvalues)
  rawpvalues <- rawpvalues[ord.inx]
  clapvalues <- p.adjust(rawpvalues, method ="fdr")

  # p.adjust
  wil_datadf <- as.data.frame(otu[,ord.inx])


  ldares <- MASS::lda(claslbl ~ .,data = wil_datadf)
  # ldares
  ldamean <- as.data.frame(t(ldares$means))
  ldamean
  class_no <<- length(unique(claslbl))
  ldamean$max <- apply(ldamean[,1:class_no],1,max);
  ldamean$min <- apply(ldamean[,1:class_no],1,min);
  #---计算LDA
  ldamean$LDAscore <- signif(log10(1+abs(ldamean$max-ldamean$min)/2),digits=3);
  head(ldamean)

  a = rep("A",length(ldamean$max))
  for (i in 1:length(ldamean$max)) {
    name =colnames(ldamean[,1:class_no])
    a[i] = name[ldamean[,1:class_no][i,] %in% ldamean$max[i]]
  }
  ldamean$class = a

  tem1 = row.names(ldamean)
  tem1 %>% as.character()
  ldamean$Pvalues <- signif(rawpvalues[match(row.names(ldamean),names(rawpvalues))],digits=5)
  ldamean$FDR <- signif(clapvalues,digits=5)
  resTable <- ldamean
  rawNms <- rownames(resTable);
  rownames(resTable) <- gsub("`", '', rawNms);


  if (adjust.p) {
    de.Num <- sum(clapvalues <= p.lvl & ldamean$LDAscore>=lda.lvl)

  } else {
    de.Num <- sum(rawpvalues <= p.lvl & ldamean$LDAscore>=lda.lvl)
  }

  if(de.Num == 0){
    current.msg <<- "No significant features were identified with given criteria.";
  }else{
    current.msg <<- paste("A total of", de.Num, "significant features with given criteria.")
  }
  print(current.msg)
  # sort by p value
  ord.inx <- order(resTable$Pvalues, resTable$LDAscore)
  resTable <- resTable[ord.inx, ,drop=FALSE]
  resTable <- resTable[,c(ncol(resTable),1:(ncol(resTable)-1))]
  resTable <- resTable[,c(ncol(resTable),1:(ncol(resTable)-1))]

  # resTable %>% tail()
  ldamean$Pvalues[is.na(ldamean$Pvalues)] = 1
  if (adjust.p) {
    taxtree = resTable[clapvalues <=p.lvl & ldamean$LDAscore>=lda.lvl,]
  } else {
    # taxtree = resTable[ldamean$Pvalues <=p.lvl & ldamean$LDAscore>=lda.lvl,]
    taxtree = resTable[ldamean$Pvalues <=p.lvl,]
  }

  #-提取所需要的颜色
  colour = c('darkgreen','red',"blue","#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
  selececol = colour[1:length(levels(as.factor(taxtree$class)))]
  names(selececol) = levels(as.factor(taxtree$class))
  A = rep("a",length(row.names(taxtree)))

  for (i in 1:length(row.names(taxtree))) {
    A[i] = selececol [taxtree$class[i]]
  }

  taxtree$color = A
  # taxtree <- taxtree[row.names(taxtree) != "k__Bacteria",]
  # node_ids <- p0$data
  # anno <- rep("white", nrow(p1$data))

  lefse_lists = data.frame(node=row.names(taxtree),
                           color=A,
                           Group = taxtree$class,
                           stringsAsFactors = FALSE
  )


  return(list(lefse_lists,taxtree))
}
