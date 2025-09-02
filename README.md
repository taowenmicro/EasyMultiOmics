# EasyMultiOmics: A Multi-Omics Analysis Framework Based on Cross-Omics Interaction Algorithm

## Main features:

*   Provides a comprehensive multi-omics analysis pipeline with emphasis on count data and mass spectrometry data integration.
*   Enables bidirectional correlation analysis between count data and mass spectrometry data.
*   Incorporates intelligent feature selection with statistical power consideration.
*   Offers sophisticated data alignment for batch effect correction and multi-omics integration.
*   Generates structured output with publication-ready visualizations and organized result directories.

## R Package Structure

*   **/data**: This folder contains example datasets used for demonstrating the functionality of the package.
*   **/R**: This folder contains the source code for all the functions in the package. Each file corresponds to a specific analysis or utility function.
*   **/pipeline**: This folder contains the workflow scripts for running different analyses. These scripts provide a structured approach to performing multi-omics data mining.

## Installation

To install EasyMultiOmics, you can use the following commands in R:

```R
# Install devtools if not already installed
install.packages("devtools")
library(devtools)

# Install EasyMultiOmics from GitHub
devtools::install_github("taowenmicro/EasyMultiOmics")
```

## Functionality Overview

### Detailed Functionality of `across_omics.R`

The `across_omics.R` script is a key component of the R package, designed to perform integrated analysis across different omics datasets. Here’s a detailed breakdown of its functionality:

#### Load Libraries

Start by loading the necessary R packages.

```R
# Load libraries
library(sva)          # ComBat batch correction
library(RSpectra)     # Sparse PCA
library(dplyr)        # Data manipulation
library(compositions) # CLR
library(ggClusterNet)
library(phyloseq)
library(devtools)
load_all()
library(EasyMultiOmics)
library(caret)
library("DESeq2")
library(limma)
library(randomForest)
library(DESeq2)
library(igraph)
library(openxlsx)
```

#### Create Main Directory Structure

Create the main directory structure for storing results.

```R
# Create main directory ----
main_dir <- "./result/across_omics"
if (!dir.exists(main_dir)) dir.create(main_dir)
```

#### Feature Selection in Genomics (Microbiome Example)

Perform feature selection in genomics data using various methods and thresholds.

```R
# Feature selection in genomics (using microbiome data as an example)#------
ps = ps.16s %>% subset_samples.wt("Group",c("WT","KO"))
# ps = ps %>% subset_taxa.wt("Species",c("unclassified"),TRUE)
ps = ps %>% tax_glom_wt(6)
ps

res = count_selection(ps,
                      group_var = "Group",
                      prevalence_threshold = 0.1,
                      detection_threshold = 0.001,
                      diff_method = "DESeq2",
                      ml_method = "rf",
                      cor_method = "spearman",
                      cor_threshold = 0.6,
                      p_threshold = 0.05,
                      weights = list(abundance = 0.3,
                                     importance = 0.3,
                                     differential = 0.2,
                                     network = 0.2))

dat = res$final_ranking %>% head(100)

# Save feature selection results
genomics_wb <- createWorkbook()
addWorksheet(genomics_wb, "genomics_features")
addWorksheet(genomics_wb, "genomics_selection_stats")
writeData(genomics_wb, "genomics_features", dat, rowNames = TRUE)
if(!is.null(res$stats)) {
  writeData(genomics_wb, "genomics_selection_stats", res$stats, rowNames = TRUE)
}
saveWorkbook(genomics_wb, file.path(main_dir, "genomics_feature_selection_results.xlsx"), overwrite = TRUE)

```

#### Finding Associated Metabolites with Target Microbe

Identify metabolites associated with a target microbe.

```R
##Find Associated Metabolites with Target Microbe------
ps2 = ps.ms %>% subset_samples.wt("Group",c("WT","OE"))
tab = ps2 %>% vegan_otu() %>%
  as.data.frame()
ps = ps.16s %>% subset_samples.wt("Group",c("WT","OE")) %>% tax_glom_wt(6)

results <- count_to_mass(
  phyloseq_obj = ps,
  metabolome_mat = tab,
  target_microbe_name = "Actinocorallia" ,
  top_n = 30
)

names(results)
results$top_metabolites
tax = ps2 %>% tax_table() %>% as.data.frame()
head(tax)
dat = results$top_metabolites %>% left_join(tax,by = c("metabolite"="metab_id"))

# Save microbe-metabolite association results
microbe_metabolite_wb <- createWorkbook()
addWorksheet(microbe_metabolite_wb, "microbe_metabolite_association")
addWorksheet(microbe_metabolite_wb, "association_details")
writeData(microbe_metabolite_wb, "microbe_metabolite_association", dat, rowNames = TRUE)
if(!is.null(results$correlation_matrix)) {
  writeData(microbe_metabolite_wb, "association_details", results$correlation_matrix, rowNames = TRUE)
}
saveWorkbook(microbe_metabolite_wb, file.path(main_dir, "microbe_metabolite_association_results.xlsx"), overwrite = TRUE)
```

#### Multi-Omics Alignment

Perform multi-omics alignment to correct and align different omics datasets.

```R
# Multi-omics alignment
otu = ps.16s %>%
  tax_glom_wt(6) %>%
  vegan_otu() %>%
  as.data.frame()

tab = ps.ms %>% vegan_otu() %>%
  as.data.frame()

omics_list = list(micro = otu,ms = tab)

alignment_results <- multi_omics_alignment(
  microbiome_data = otu,  # Sample x Microbiome
  metabolome_data = tab,  # Sample x Microbiome
  n_factors = 10,                       # Number of latent factors
  method = "fast"                       # Automatically choose between fast and complete methods
)

# Save alignment results
alignment_wb <- createWorkbook()
addWorksheet(alignment_wb, "alignment_quality")
addWorksheet(alignment_wb, "aligned_microbiome")
addWorksheet(alignment_wb, "aligned_metabolome")
writeData(alignment_wb, "alignment_quality", alignment_results$quality_details, rowNames = TRUE)
writeData(alignment_wb, "aligned_microbiome", alignment_results$aligned_microbiome, rowNames = TRUE)
writeData(alignment_wb, "aligned_metabolome", alignment_results$aligned_metabolome, rowNames = TRUE)
saveWorkbook(alignment_wb, file.path(main_dir, "multi_omics_alignment_results.xlsx"), overwrite = TRUE)
```

### Other Pipeline Scripts (Example: `1.1.pipeline.metm.R`)

Other pipeline scripts in the package follow a similar structure, focusing on specific types of analyses. For example, the `1.1.pipeline.metm.R` script performs microbiome analysis.

#### Read Data and Set Parameters

Start by reading the necessary data and setting up parameters.

```R
rm(list=ls())

# BiocManager::install("MicrobiotaProcess")
library(EasyMultiOmics)
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(ggrepel)
library(openxlsx)

# Effective Gene Count
sample_sums(ps.micro)

## Set parameters-----
map= sample_data(ps.micro)
head(map)
phyloseq::tax_table(ps.micro) %>% head()

# Extract the number of unique group factors
gnum = phyloseq::sample_data(ps.micro)$Group %>% unique() %>% length()
gnum

#--Set the order for plotting according to the order in the map file of the ps.micro object
axis_order =  phyloseq::sample_data(ps.micro)$Group %>%unique();axis_order
col.g =  c("KO" = "#D55E00", "WT" = "#0072B2", "OE" = "#009E73")

#-theme--
package.amp()
res = theme_my(ps.micro)
mytheme1 = res[[1]]
mytheme2 = res[[2]];
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]
```

#### Data Analysis

Perform various analyses, such as alpha and beta diversity calculations.

```R
# alpha diversity -----
#1 alpha.metm:#----
all.alpha = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )

#--Calculate alpha diversity metrics
tab = alpha.metm(ps = ps.micro,group = "Group" )
head(tab)
data = cbind(data.frame(ID = 1:length(tab$Group),group = tab$Group),tab[all.alpha])
head(data)
data$ID = as.character(data$ID)

result = MuiKwWlx2(data = data,num = 3:6)

result1 = FacetMuiPlotresultBox(data = data,num = 3:6,
                                          result = result,
                                          sig_show ="abc",ncol = 4,width = 0.4 )


p1_1 = result1[[1]] +scale_fill_manual(values = col.g)

p1_1+
  ggplot2::scale_x_discrete(limits = axis_order) +
  # theme_cell()+
  theme_nature()+

  ggplot2::guides(fill = guide_legend(title = none))

# Save path
alppath <- file.path(path, "alpha")
dir.create(alppath, showWarnings = FALSE)

# Save Alpha diversity images
ggsave(file.path(alppath, "alpha_diversity_boxplot.png"), plot = p1_1, width = 10, height = 8)
ggsave(file.path(alppath, "alpha_diversity_boxplot.pdf"), plot = p1_1, width = 10, height = 8)
```

