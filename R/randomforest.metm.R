
#' @title Random forest modeling for microbiome data
#' @description Random forest modeling for microbiome data
#' @param otu OTU/ASV table;
#' @param map Sample metadata;
#' @param tax taxonomy table;
#' @param ps phyloseq object of microbiome;
#' @param Group column name for groupID in map table.
#' @param optimal important OTU number which selected
#' @param rfcv TURE or FELSE,whether need to do cross-validation
#' @param nrfcvnum Number of cross-validation
#' @param min Circle diagram inner diameter adjustment
#' @param max Circle diagram outer diameter adjustment
#' @return list contain ggplot object and table.
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' result = randomforest.metm(ps = ps,group  = "Group",optimal = 20,rfcv = TRUE,nrfcvnum = 5,min = -1,max = 5)
#'@export

randomforest.metm =function (otu = NULL, tax = NULL, map = NULL, tree = NULL, ps = NULL,
          group = "Group", optimal = 20, rfcv = FALSE, nrfcvnum = 5,
          min = -1, max = 5, fill = "Phylum", lab = "id")
{
  ps = ggClusterNet::inputMicro(otu, tax, map, tree, ps, group = group) %>%
    ggClusterNet::scale_micro()
  map = as.data.frame(phyloseq::sample_data(ps))
  mapping = as.data.frame(phyloseq::sample_data(ps))
  otutab = as.data.frame((ggClusterNet::vegan_otu(ps)))
  colnames(otutab) = paste("wentao", colnames(otutab), sep = "")
  otutab$group = factor(mapping$Group)
  colnames(otutab) <- gsub("-", "_", colnames(otutab))
  colnames(otutab) <- gsub("[/]", "_", colnames(otutab))
  colnames(otutab) <- gsub("[(]", "_", colnames(otutab))
  colnames(otutab) <- gsub("[)]", "_", colnames(otutab))
  colnames(otutab) <- gsub("[:]", "_", colnames(otutab))
  colnames(otutab) <- gsub("[[]", "_", colnames(otutab))
  colnames(otutab) <- gsub("[]]", "_", colnames(otutab))
  colnames(otutab) <- gsub("[#]", "_", colnames(otutab))
  colnames(otutab) <- gsub("[+]", "_", colnames(otutab))
  colnames(otutab) <- gsub(" ", "_", colnames(otutab))
  colnames(otutab) <- gsub("[,]", "_", colnames(otutab))
  model_rf = randomForest::randomForest(group ~ ., data = otutab,
                                        importance = TRUE, proximity = TRUE)
  print(model_rf)
  Confusion_matrix <- as.data.frame(model_rf$confusion)
  Confusion_matrix$class.error <- round(Confusion_matrix$class.error,
                                        3)
  Confusion_matrix$Group = row.names(Confusion_matrix)
  Confusion_matrix <- dplyr::select(Confusion_matrix, Group,
                                    everything())
  model_Accuracy_rates <- paste(round(100 - tail(model_rf$err.rate[,
                                                                   1], 1) * 100, 2), "%", sep = "")
  model_Accuracy_rates = data.frame(ID = "model Accuracy rates",
                                    model_Accuracy_rates = model_Accuracy_rates)
  colnames(model_Accuracy_rates) = c("Random foreest", "Fu wilt model")
  tab2 <- ggpubr::ggtexttable(Confusion_matrix, rows = NULL)
  tab1 <- ggpubr::ggtexttable(model_Accuracy_rates, rows = NULL)
  library(patchwork)
  pn <- tab1/tab2
  a = as.data.frame(round(randomForest::importance(model_rf),
                          2))
  a$id = row.names(a)
  row.names(a) = gsub("wentao", "", row.names(a))
  a$id = gsub("wentao", "", a$id)
  a2 <- dplyr::arrange(a, desc(MeanDecreaseAccuracy)) %>% as.data.frame()
  row.names(a2) = a2$id
  a3 = head(a2, n = optimal)
  head(a3)
  if (is.null(ps@tax_table)) {
    tax = data.frame(row.names = gsub("wentao", "", colnames(otutab)[-(length(colnames(otutab)))]),
                     ID = gsub("wentao", "", colnames(otutab)[-(length(colnames(otutab)))]),
                     class = gsub("wentao", "", colnames(otutab)[-(length(colnames(otutab)))]))
  }
  else {
    tax <- as.data.frame(ggClusterNet::vegan_tax(ps))
  }
  tax$org.id = row.names(tax)
  row.names(tax) = gsub("wentao", "", colnames(otutab)[-(length(colnames(otutab)))])
  tax = tax[rownames(a3), ]
  dim(a3)
  tax$id = NULL
  a3 = merge(a3, tax, by = "row.names", all = F)
  row.names(a3) = a3$Row.names
  a3$Row.names = NULL
  if (is.null(ps@tax_table)) {
    tax = data.frame(row.names = gsub("wentao", "", colnames(otutab)[-(length(colnames(otutab)))]),
                     ID = gsub("wentao", "", colnames(otutab)[-(length(colnames(otutab)))]),
                     class = gsub("wentao", "", colnames(otutab)[-(length(colnames(otutab)))]))
  }
  else {
    tax <- as.data.frame(ggClusterNet::vegan_tax(ps))
  }
  tax$org.id = row.names(tax)
  row.names(tax) = gsub("wentao", "", colnames(otutab)[-(length(colnames(otutab)))])
  tax = tax[rownames(a2), ]
  dim(a2)
  tax$id = NULL
  a2 = merge(a2, tax, by = "row.names", all = F)
  row.names(a2) = a3$Row.names
  a2$Row.names = NULL
  OTU = ggClusterNet::vegan_otu(ps)
  design = as.data.frame(phyloseq::sample_data(ps))
  iris.split <- split(as.data.frame(OTU), as.factor(design$Group))
  iris.apply <- lapply(iris.split, function(x) colMeans(x,
                                                        na.rm = TRUE))
  norm2 <- do.call(rbind, iris.apply) %>% t()
  colnames(norm2) = paste(colnames(norm2), "mean", sep = "")
  ind_fal = merge(a3, norm2, by = "row.names", all = F)
  head(ind_fal)
  if (is.null(lab)) {
    lab = "id"
  }
  if (is.null(fill)) {
    a3$fill = "class"
    fill = "fill"
  }
  tem = a3 %>% arrange(desc(MeanDecreaseAccuracy))
  head(tem)
  a3$id = factor(a3$id, level = tem$id[length(tem$id):1])
  p1 <- ggplot(a3, aes(x = MeanDecreaseAccuracy, y = id, fill = !!sym(fill),
                       color = !!sym(fill))) + geom_point(size = 6, pch = 21) +
    geom_segment(aes(yend = id), xend = 0, size = 3) + geom_text(aes(x = MeanDecreaseAccuracy *
                                                                       1.1, label = !!sym(lab)), size = 3, color = "black")
  a3 <- dplyr::arrange(a3, desc(MeanDecreaseAccuracy))
  a3$iid = paste(1:length(a3$id))
  angle1 = 90 - 360 * (as.numeric(a3$iid) - 0.5)/length(a3$id)
  a3$id = factor(a3$id, levels = a3$id)
  p2 = a3 %>% ggplot(aes(x = factor(id), y = MeanDecreaseAccuracy,
                         label = id)) + geom_bar(stat = "identity", position = "dodge",
                                                 fill = "blue") + geom_text(hjust = 0, angle = angle1,
                                                                            alpha = 1) + coord_polar() + ggtitle("") + ylim(c(min,
                                                                                                                              max)) + theme_void()
  p2
  return(list(p1, p2, a3, pn, a2))
}
