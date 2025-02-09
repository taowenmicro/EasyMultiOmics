#' @title  Perform differential analysis using ANCOM
#' @description This function is used to analyze transcriptome functional composition
#'  data  by conducting ANCOM analysis to detect significant differences in genes
#'  abundance between different groups.
#' @param ps A phyloseq format file used as an alternative for the input containing transcriptome functional composition table, tax, and sample metadata.
#' @param group Column name for groupID in map table(sample metadata).
#' @param p_adj_method Select the method for p-value correction("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#' "fdr"),when padjust is `TRUE`. Default is `BH`.
#' @param ANCOM.value Value used for ANCOM analysis, default is "detected_0.6"
#' @param alpha Significance level, default is 0.05.
#' @return A list containing the results of significant differences in gene abundance between different groups.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' \dontrun{
#' data(ps.trans)
#' ps = ps.trans%>% filter_taxa(function(x) sum(x ) > 5 , TRUE)
#' dat <- ancom.trans(ps = ps%>% filter_OTU_ps(50), group = "Group",  p_adj_method = "BH",ANCOM.value = "detected_0.6",alpha=0.05)
#' dat1 = dat$WT_OE
#' head(dat1)
#' dat2 = dat$WT_KO
#' head(dat2)
#' dat3 = dat$OE_KO
#' head(dat3)
#' }
#' @export


ancom.trans = function(
    ps = ps0,
    group = "Group",
    p_adj_method = "BH",
    ANCOM.value = "detected_0.6",
    alpha=0.05
){

  map = sample_data(ps)
  head(map)
  id.g = map$Group %>% unique() %>% as.character() %>% combn(2)

  aaa = id.g
  data_list =  list()
  i=3
  for (i in 1:dim(aaa)[2]) {
    Desep_group = aaa[,i]
    print( Desep_group)
    ps.cs = ps %>% subset_samples.wt("Group" ,id.g[,i])

    ps.cs
    ASV_table = ps.cs %>%
      # filter_taxa(function(x) sum(x ) > 200 , TRUE) %>%
      # subset_samples(Group %in% id.g[,i]) %>%
      vegan_otu() %>%round() %>% t() %>%
      as.data.frame()
    groupings <- ps.cs %>%
      # subset_samples(Group%in% id.g[,i]) %>%
      sample_data()
    groupings$ID = row.names(groupings)
    groupings$Sample <- rownames(groupings)

    prepro <- feature_table_pre_process(feature_table = ASV_table,
                                        meta_data = groupings,
                                        sample_var = 'Sample',
                                        group_var = NULL,
                                        out_cut = 0.05,
                                        zero_cut = 0.90,
                                        lib_cut = 1000,
                                        neg_lb=FALSE)

    feature_table <- prepro$feature_table
    metadata <- prepro$meta_data
    struc_zero <- prepro$structure_zeros

    #run ancom

    adj_formula=NULL
    rand_formula=NULL
    res <- ANCOM(feature_table = feature_table,
                 meta_data = metadata,
                 struc_zero = struc_zero,
                 main_var = group,
                 p_adj_method = p_adj_method,
                 alpha=alpha, adj_formula = adj_formula,
                 rand_formula = rand_formula
    )
    dat = res$out
    head(dat)
    tab.d2 = dat %>% dplyr::select(taxa_id,detected_0.6) %>%
      dplyr::filter(detected_0.6 == TRUE) %>%
      dplyr::rename(
        OTU = taxa_id,
        p = detected_0.6
      )  %>%
      dplyr::mutate(group = "ANCOMII")

    # head(tab.d2)
    # 不同阈值的显著差异的数量不同，随着阈值越来越大，显著的OTU也越来越少
    # write.table(res$out, file= "./ANCOM.txt", quote=FALSE, sep="\t", col.names = NA)
    data_list[[i]] = data.frame(micro = tab.d2$OTU,method = tab.d2$group,adjust.p = tab.d2$p)

    names( data_list)[i]= Desep_group %>% paste( collapse = "_")
  }
  return(data_list)
}


feature_table_pre_process = function(feature_table, meta_data, sample_var, group_var = NULL,
                                     out_cut = 0.05, zero_cut = 0.90, lib_cut, neg_lb){
  feature_table = data.frame(feature_table, check.names = FALSE)
  meta_data = data.frame(meta_data, check.names = FALSE)
  # Drop unused levels
  meta_data[] = lapply(meta_data, function(x) if(is.factor(x)) factor(x) else x)
  # Match sample IDs between metadata and feature table
  sample_ID = intersect(meta_data[, sample_var], colnames(feature_table))
  feature_table = feature_table[, sample_ID]
  meta_data = meta_data[match(sample_ID, meta_data[, sample_var]), ]

  # 1. Identify outliers within each taxon
  if (!is.null(group_var)) {
    group = meta_data[, group_var]
    z = feature_table + 1 # Add pseudo-count (1)
    f = log(z); f[f == 0] = NA; f = colMeans(f, na.rm = T)
    f_fit = lm(f ~ group)
    e = rep(0, length(f)); e[!is.na(group)] = residuals(f_fit)
    y = t(t(z) - e)

    outlier_check = function(x){
      # Fitting the mixture model using the algorithm of Peddada, S. Das, and JT Gene Hwang (2002)
      mu1 = quantile(x, 0.25, na.rm = T); mu2 = quantile(x, 0.75, na.rm = TRUE)
      sigma1 = quantile(x, 0.75, na.rm = T) - quantile(x, 0.25, na.rm = TRUE); sigma2 = sigma1
      pi = 0.75
      n = length(x)
      epsilon = 100
      tol = 1e-5
      score = pi*dnorm(x, mean = mu1, sd = sigma1)/((1 - pi)*dnorm(x, mean = mu2, sd = sigma2))
      while (epsilon > tol) {
        grp1_ind = (score >= 1)
        mu1_new = mean(x[grp1_ind]); mu2_new = mean(x[!grp1_ind])
        sigma1_new = sd(x[grp1_ind]); if(is.na(sigma1_new)) sigma1_new = 0
        sigma2_new = sd(x[!grp1_ind]); if(is.na(sigma2_new)) sigma2_new = 0
        pi_new = sum(grp1_ind)/n

        para = c(mu1_new, mu2_new, sigma1_new, sigma2_new, pi_new)
        if(any(is.na(para))) break

        score = pi_new * dnorm(x, mean = mu1_new, sd = sigma1_new)/
          ((1-pi_new) * dnorm(x, mean = mu2_new, sd = sigma2_new))

        epsilon = sqrt((mu1 - mu1_new)^2 + (mu2 - mu2_new)^2 +
                         (sigma1 - sigma1_new)^2 + (sigma2 - sigma2_new)^2 + (pi - pi_new)^2)
        mu1 = mu1_new; mu2 = mu2_new; sigma1 = sigma1_new; sigma2 = sigma2_new; pi = pi_new
      }

      if(mu1 + 1.96 * sigma1 < mu2 - 1.96 * sigma2){
        if(pi < out_cut){
          out_ind = grp1_ind
        }else if(pi > 1 - out_cut){
          out_ind = (!grp1_ind)
        }else{
          out_ind = rep(FALSE, n)
        }
      }else{
        out_ind = rep(FALSE, n)
      }
      return(out_ind)
    }
    out_ind = matrix(FALSE, nrow = nrow(feature_table), ncol = ncol(feature_table))
    out_ind[, !is.na(group)] = t(apply(y, 1, function(i)
      unlist(tapply(i, group, function(j) outlier_check(j)))))

    feature_table[out_ind] = NA
  }

  # 2. Discard taxa with zeros  >=  zero_cut
  zero_prop = apply(feature_table, 1, function(x) sum(x == 0, na.rm = TRUE)/length(x[!is.na(x)]))
  taxa_del = which(zero_prop >= zero_cut)
  if(length(taxa_del) > 0){
    feature_table = feature_table[- taxa_del, ]
  }

  # 3. Discard samples with library size < lib_cut
  lib_size = colSums(feature_table, na.rm = TRUE)
  if(any(lib_size < lib_cut)){
    subj_del = which(lib_size < lib_cut)
    feature_table = feature_table[, - subj_del]
    meta_data = meta_data[- subj_del, ]
  }

  # 4. Identify taxa with structure zeros
  if (!is.null(group_var)) {
    group = factor(meta_data[, group_var])
    present_table = as.matrix(feature_table)
    present_table[is.na(present_table)] = 0
    present_table[present_table != 0] = 1

    p_hat = t(apply(present_table, 1, function(x)
      unlist(tapply(x, group, function(y) mean(y, na.rm = T)))))
    samp_size = t(apply(feature_table, 1, function(x)
      unlist(tapply(x, group, function(y) length(y[!is.na(y)])))))
    p_hat_lo = p_hat - 1.96 * sqrt(p_hat * (1 - p_hat)/samp_size)

    struc_zero = (p_hat == 0) * 1
    # Whether we need to classify a taxon into structural zero by its negative lower bound?
    if(neg_lb) struc_zero[p_hat_lo <= 0] = 1

    # Entries considered to be structural zeros are set to be 0s
    struc_ind = struc_zero[, group]
    feature_table = feature_table * (1 - struc_ind)

    colnames(struc_zero) = paste0("structural_zero (", colnames(struc_zero), ")")
  }else{
    struc_zero = NULL
  }

  # 5. Return results
  res = list(feature_table = feature_table, meta_data = meta_data, structure_zeros = struc_zero)
  return(res)
}

# ANCOM main function
ANCOM = function(feature_table, meta_data, struc_zero = NULL, main_var, p_adj_method = "BH",
                 alpha = 0.05, adj_formula = NULL, rand_formula = NULL, ...){
  # OTU table transformation:
  # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
  if (!is.null(struc_zero)) {
    num_struc_zero = apply(struc_zero, 1, sum)
    comp_table = feature_table[num_struc_zero == 0, ]
  }else{
    comp_table = feature_table
  }
  comp_table = log(as.matrix(comp_table) + 1)
  n_taxa = dim(comp_table)[1]
  taxa_id = rownames(comp_table)
  n_samp = dim(comp_table)[2]

  # Determine the type of statistical test and its formula.
  if (is.null(rand_formula) & is.null(adj_formula)) {
    # Basic model
    # Whether the main variable of interest has two levels or more?
    if (length(unique(meta_data%>%pull(main_var))) == 2) {
      # Two levels: Wilcoxon rank-sum test
      tfun = exactRankTests::wilcox.exact
    } else{
      # More than two levels: Kruskal-Wallis test
      tfun = stats::kruskal.test
    }
    # Formula
    tformula = formula(paste("x ~", main_var, sep = " "))
  }else if (is.null(rand_formula) & !is.null(adj_formula)) {
    # Model: ANOVA
    tfun = stats::aov
    # Formula
    tformula = formula(paste("x ~", main_var, "+", adj_formula, sep = " "))
  }else if (!is.null(rand_formula)) {
    # Model: Mixed-effects model
    tfun = nlme::lme
    # Formula
    if (is.null(adj_formula)) {
      # Random intercept model
      tformula = formula(paste("x ~", main_var))
    }else {
      # Random coefficients/slope model
      tformula = formula(paste("x ~", main_var, "+", adj_formula))
    }
  }

  # Calculate the p-value for each pairwise comparison of taxa.
  p_data = matrix(NA, nrow = n_taxa, ncol = n_taxa)
  colnames(p_data) = taxa_id
  rownames(p_data) = taxa_id
  for (i in 1:(n_taxa - 1)) {
    # Loop through each taxon.
    # For each taxon i, additive log ratio (alr) transform the OTU table using taxon i as the reference.
    # e.g. the first alr matrix will be the log abundance data (comp_table) recursively subtracted
    # by the log abundance of 1st taxon (1st column) column-wisely, and remove the first i columns since:
    # the first (i - 1) columns were calculated by previous iterations, and
    # the i^th column contains all zeros.
    alr_data = apply(comp_table, 1, function(x) x - comp_table[i, ])
    # apply(...) allows crossing the data in a number of ways and avoid explicit use of loop constructs.
    # Here, we basically want to iteratively subtract each column of the comp_table by its i^th column.
    alr_data = alr_data[, - (1:i), drop = FALSE]
    n_lr = dim(alr_data)[2] # number of log-ratios (lr)
    alr_data = cbind(alr_data, meta_data) # merge with the metadata

    # P-values
    if (is.null(rand_formula) & is.null(adj_formula)) {
      p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
        tfun(tformula, data = data.frame(x, alr_data, check.names = FALSE))$p.value
      }
      )
    }else if (is.null(rand_formula) & !is.null(adj_formula)) {
      p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
        fit = tfun(tformula,
                   data = data.frame(x, alr_data, check.names = FALSE),
                   na.action = na.omit)
        summary(fit)[[1]][main_var, "Pr(>F)"]
      }
      )
    }else if (!is.null(rand_formula)) {
      p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
        fit = tfun(fixed = tformula,
                   data = data.frame(x, alr_data, check.names = FALSE),
                   random = formula(rand_formula),
                   na.action = na.omit, ...)
        anova(fit)[main_var, "p-value"]
      }
      )
    }
  }
  # Complete the p-value matrix.
  # What we got from above iterations is a lower triangle matrix of p-values.
  p_data[upper.tri(p_data)] = t(p_data)[upper.tri(p_data)]
  diag(p_data) = 1 # let p-values on diagonal equal to 1

  # Multiple comparisons correction.
  q_data = apply(p_data, 2, function(x) p.adjust(x, method = p_adj_method))

  # Calculate the W statistic of ANCOM.
  # For each taxon, count the number of q-values < alpha.
  W = apply(q_data, 2, function(x) sum(x < alpha))

  # Organize outputs
  out_comp = data.frame(taxa_id, W, row.names = NULL, check.names = FALSE)
  # Declare a taxon to be differentially abundant based on the quantile of W statistic.
  # We perform (n_taxa - 1) hypothesis testings on each taxon, so the maximum number of rejections is (n_taxa - 1).
  out_comp = out_comp%>%mutate(detected_0.9 = ifelse(W > 0.9 * (n_taxa -1), TRUE, FALSE),
                               detected_0.8 = ifelse(W > 0.8 * (n_taxa -1), TRUE, FALSE),
                               detected_0.7 = ifelse(W > 0.7 * (n_taxa -1), TRUE, FALSE),
                               detected_0.6 = ifelse(W > 0.6 * (n_taxa -1), TRUE, FALSE))

  # Taxa with structural zeros are automatically declared to be differentially abundant
  if (!is.null(struc_zero)){
    out = data.frame(taxa_id = rownames(struc_zero), W = Inf, detected_0.9 = TRUE,
                     detected_0.8 = TRUE, detected_0.7 = TRUE, detected_0.6 = TRUE,
                     row.names = NULL, check.names = FALSE)
    out[match(taxa_id, out$taxa_id), ] = out_comp
  }else{
    out = out_comp
  }

  # library(contrib)
  # Draw volcano plot
  # Calculate clr
  # library(edgeR)
  # clr_table = apply(feature_table, 2, clr)
  clr_table = feature_table
  # Calculate clr mean difference
  eff_size = apply(clr_table, 1, function(y)
    lm(y ~ x, data = data.frame(y = y,
                                x = meta_data %>% pull(main_var),
                                check.names = FALSE))$coef[-1])

  if (is.matrix(eff_size)){
    # Data frame for the figure
    dat_fig = data.frame(taxa_id = out$taxa_id, t(eff_size), y = out$W, check.names = FALSE) %>%
      mutate(zero_ind = factor(ifelse(is.infinite(y), "Yes", "No"), levels = c("Yes", "No"))) %>%
      gather(key = group, value = x, rownames(eff_size))
    # Replcace "x" to the name of covariate
    dat_fig$group = sapply(dat_fig$group, function(x) gsub("x", paste0(main_var, " = "), x))
    # Replace Inf by (n_taxa - 1) for structural zeros
    dat_fig$y = replace(dat_fig$y, is.infinite(dat_fig$y), n_taxa - 1)

    fig = ggplot(data = dat_fig) + aes(x = x, y = y) +
      geom_point(aes(color = zero_ind)) +
      facet_wrap(~ group) +
      labs(x = "CLR mean difference", y = "W statistic") +
      scale_color_discrete(name = "Structural zero", drop = FALSE) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "top",
            strip.background = element_rect(fill = "white"))
    fig
  } else{
    # Data frame for the figure
    dat_fig = data.frame(taxa_id = out$taxa_id, x = eff_size, y = out$W) %>%
      mutate(zero_ind = factor(ifelse(is.infinite(y), "Yes", "No"), levels = c("Yes", "No")))
    # Replace Inf by (n_taxa - 1) for structural zeros
    dat_fig$y = replace(dat_fig$y, is.infinite(dat_fig$y), n_taxa - 1)

    fig = ggplot(data = dat_fig) + aes(x = x, y = y) +
      geom_point(aes(color = zero_ind)) +
      labs(x = "CLR mean difference", y = "W statistic") +
      scale_color_discrete(name = "Structural zero", drop = FALSE) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "top")
    fig
  }

  res = list(out = out, fig = fig)
  return(res)
}






