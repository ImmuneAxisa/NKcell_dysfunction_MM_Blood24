#' title: "Hitchhiker's guide to NK"

##########################
##   Loading packages   ##
##########################
suppressPackageStartupMessages({
  library(tidyverse)
  library(conflicted)
  library(qs)
  library(Seurat)
  library(harmony)
  library(clustree)
})
################################################################################
## Starting with the loading of the end product of "04_QC_check_and_trim.R"   ##
################################################################################

trim_sc <- qread("data/qc_trimmed_t_plus_nk.qs",nthreads = 10)

preprocess_plots<- "../pub_ready/docs/preprocess_plots/"
dir.create(path = preprocess_plots)
#### Preprocessing within Seurat ####
trim_sc <- NormalizeData(trim_sc,
                     normalization.method = "LogNormalize",
                     scale.factor = 10000) %>% FindVariableFeatures() %>% ScaleData(vars.to.regress=c("percent.mt","percent.rb")) %>% RunPCA()

#' After normalization, finding most variable features, scaling for those, and then PCs of the dataset is calculated.


trim_sc <- RunHarmony(trim_sc,reduction = "pca",group.by.vars = "id",verbose = T)
# Harmony package is utilized with all 50 PC dims.
# To integrate different scRNAseq emulsion runs of this dataset and reducing the technical variance between samples while keeping biological variance of the samples
#'
# **Selection of optimum number of dimensions**
source(file = "../pub_ready/func/quant_harmonys.R")
#' The helper function above checks two criteria; how many dimensions are required to cover 90% of the cumulative standard deviation of the dataset,
#' and when the difference of standard deviation between two consecutive dimensions drops less than 0.1 percent.
#' Then ranking the minimum dimension of those two, to continue analysis as the optimum dimension to use downstream.

#'
#' Below is the ranking of the harmony reduced dimension of the dataset and the selection of optimum dims of harmony reduction for clustering and UMAP.
p <- quant_opt_harmonys(trim_sc)
# Saving this plot
ggsave(
  filename = "optimum_number_of_harmony_dims.png",
  p,
  device = "png",
  dpi = 150,
  path = preprocess_plots,
  height = 8,
  width = 10)

# According to the quantitative assessment of the harmony dimensional reductions of the dataset,
#  using up to **30th harmony** on downstream calculations is appropriate.
trim_sc <- FindNeighbors(trim_sc,reduction = "harmony",dims = 1:30,k.param = 20)

# After determining the neighbors of the cells,
# we are going to determine unsupervised clusters of cells within a neighborhood for a set of sequential granularity.
# This will help us to better appreciate the migratory routes cells within clusters for different resolutions and distinct origins of clusters.

## Clustering by Leiden algorithm by sequential resolutions ##
trim_sc <- FindClusters(trim_sc,algorithm = 4,method = "igraph",resolution = seq(0.1,1,0.1),random.seed = 2206)

# Calculating Umap Coordinates of the cells
trim_sc <- RunUMAP(trim_sc,reduction = "harmony",dims = 1:30)

# Saving up the preprocess completed dataset
qsave(trim_sc,file = "data/preprocessed_t_nk.qs",nthreads = 10)
