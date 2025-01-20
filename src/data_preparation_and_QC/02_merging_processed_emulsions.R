suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(SeuratDisk)
  library(SeuratWrappers)
  library(qs)
})

# Load all samples and merge them into a single Seurat object
name_of_samples <- list.files(path = "data/velo_embedded_wnn/")
# This "Merge" part is quite heavy in memory !
# if you have LESS than *100 GB* of RAM
# use this;
#         big_sc <- lapply(name_of_samples, function(file) {LoadH5Seurat(file = paste0("data/velo_embedded_wnn/", file),assays = c("RNA", "ADT"))})
# otherwise, ignore this comment
big_sc <- lapply(name_of_samples, function(file) {
  LoadH5Seurat(file = paste0("data/velo_embedded_wnn/", file),
               assays = c("RNA", "ADT","spliced","unspliced","ambigous"))
})
big_sc <- Reduce(merge, big_sc)
gc()

# Save the resulting Seurat object
qsave(x = big_sc, file = "velo_embedded_big_sc_whole_pbmc_annotated.qs")
