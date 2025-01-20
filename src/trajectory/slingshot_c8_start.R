library(qs)
library(Seurat)
library(slingshot)
setwd(here::here())

sc_object <- qread("results/NK_analysis/seurat_objects/analysisNK_Ruchan.qs")


rd <- Embeddings(object = sc_object, reduction = "harmony")
Clusters <- sc_object$seurat_clusters
colnames(rd) <- colnames(Embeddings(object = sc_object, reduction = "harmony"))
df <- data.frame(rd,Embeddings(object = sc_object, reduction = "umap"), "Clusters" = as.character(Clusters))



##### Density Shrink, do not allow tree break
sds <-
  slingshot(
    rd,
    clusterLabels = Clusters,
    start.clus = 8,
    shrink.method = "density",
    allow.breaks = TRUE) 



qsave(sds, "results/NK_analysis/trajectory_pseudotime/slingshot_pto_c8.qs")