library(qs)
library(magrittr)
library(Seurat)
library(slingshot)
setwd(here::here())

sc_object <- qread("results/NK_analysis/seurat_objects/analysisNK_Ruchan.qs")

slingshot_wrapper <- function(sc_data, start_clust) {
  
  rd <- Embeddings(object = sc_data, reduction = "harmony")
  Clusters <- sc_data$seurat_clusters
  colnames(rd) <- colnames(Embeddings(object = sc_data, reduction = "harmony"))
  df <- data.frame(rd,Embeddings(object = sc_data, reduction = "umap"), "Clusters" = as.character(Clusters))
  
  
  
  ##### Density Shrink, do not allow tree break
  sds <-
    slingshot(
      rd,
      clusterLabels = Clusters,
      start.clus = start_clust,
      shrink.method = "density",
      allow.breaks = TRUE) 
  
}


sc_object[, sc_object$source == "PBMC"] %>%
  slingshot_wrapper(5) %>%
  qsave("results/NK_analysis/trajectory_pseudotime/slingshot_pto_c5_pbmc.qs")

sc_object[, sc_object$source == "BMMC"] %>%
  slingshot_wrapper(5) %>%
  qsave("results/NK_analysis/trajectory_pseudotime/slingshot_pto_c5_bmmc.qs")

