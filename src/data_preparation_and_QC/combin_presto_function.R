suppressPackageStartupMessages({
 library(qs)
 library(Seurat)
 library(presto)
 library(msigdbr)
 library(fgsea)
 library(dplyr)
 library(ggplot2)
 library(tidyverse)
 library(scran)
 library(here)
 library(tictoc)
})

combin_presto <- function(sc_obj) {
 if ("curated_clusters" %in% names(sc_obj@meta.data) == FALSE)
 {sc_obj$curated_clusters <- sc_obj$RNA_snn_res.0.5
 curated_clusters <- sc_obj@meta.data$curated_clusters}
 else
 {curated_clusters <- sc_obj@meta.data$curated_clusters
 }

 combinations <-
  combn(unique(curated_clusters),
        m = 2,
        simplify = F)
 {
  tic()
  paired_results <- lapply(combinations, function(clus_comb) {
   wilcoxauc(
    sc_obj,
    "curated_clusters",
    groups_use = clus_comb,
    assay = "data"
   ) %>%
    mutate(groupA = clus_comb[1]) %>%
    mutate(groupB = clus_comb[2]) %>%
    mutate(pairs = paste0("Cl_", clus_comb[2], "vs", clus_comb[1])) %>%
    subset(group == groupB)  # Removing the doubled results of "A vs B" and "B vs A"
  })
  # Order of Cl_clus_comb[2] vs clus_comb[1] is
  # validated by checking the expected expression comparisons for clusters
  toc()
  return(paired_results)
  }
}
