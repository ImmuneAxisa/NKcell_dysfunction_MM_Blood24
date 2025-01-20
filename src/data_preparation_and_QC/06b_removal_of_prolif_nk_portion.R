##########
library(qs)
library(Seurat)
library(harmony)
library(leiden)
nk <- qread(file = "denoised_wo_velo_nk")
###
non_pro_nk <- subset(nk, idents = 1:8) # This subsetting includes Proliferating NK inferred portion of the cell,
# plus unsupervised cluster of majority of it classified as "Proliferating Nk" cells
non_pro_nk <-
 DietSeurat(
  object = non_pro_nk,
  counts = T,
  scale.data = F,
  dimreducs = NULL,
  graphs = NULL
 )

non_pro_nk <-
 non_pro_nk %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData(vars.to.regress = c("percent.mt", "percent.rb"))  %>% RunPCA()
non_pro_nk <- RunHarmony(object = non_pro_nk, group.by.vars = "id")
ElbowPlot(non_pro_nk, reduction = "harmony") # n dim 15

non_pro_nk <-
 FindNeighbors(
  non_pro_nk,
  reduction = "harmony",
  dims = 1:20,
  force.recalc = T,
  assay = "RNA"
 )
non_pro_nk <-
 FindClusters(
  object = non_pro_nk,
  algorithm = "leiden",
  method = "igraph",
  resolution = seq(0.1, 0.1, 1),
  # resolution= 0.5,
  random.seed = 2206,
  verbose = T
 )


non_pro_nk <- RunUMAP(
 non_pro_nk,
 reduction = "harmony",
 dims = 1:20,
 umap.method = "umap-learn",
 metric = "correlation",
 min.dist = 0.01,
 spread = 1
)
####
Idents(non_pro_nk) <- non_pro_nk$RNA_snn_res.0.5 # Decided to stick to this granularity level of clustering
qsave(non_pro_nk, file = "nonprolif_NKs_multiple_res_clusters_20dec22.qs")


mast_markers <-
 FindAllMarkers(
  non_pro_nk,
  assay = "RNA",
  random.seed = 2206,
  max.cells.per.ident = 1000,
  verbose = T,
  only.pos = F,
  test.use = "MAST"
 )

mast_markers <-
 mast_markers %>% scCustomize::Add_Pct_Diff(, overwrite = T) %>% left_join(,
                                                                           y = unique(annotations[, c("gene_name", "description")]),
                                                                           by = c("gene" = "gene_name"))

write.table(
 mast_markers,
 file = "Cluster_markers_by_MAST_w_percent_diffs_for_8_cluster.txt",
 quote = F,
 sep = "\t",
 row.names = F
)

ann <-
 read_rds(file = "/mnt/SERVER-CRCT-STORAGE/CRCT13/PPA/PPA_work/project/01_HCvsMM/01A_HCvsMM_10x/data/features/gene_annotations_positions.rds")
ann <-
 ann %>% dplyr::filter(gene_biotype == "protein_coding") %>% left_join(, y = unique(annotations[, c("gene_name", "description")]), by = "gene_name")
qsave(ann, file = "gene_annotations_of_this_run.qs")


nk <-
 qread(file = "data/nonprolif_NKs_multiple_res_clusters.qs")
SCpubr::do_ChordDiagramPlot(
 sample = nk,
 from = "source",
 to = "seurat_clusters",
 colors.from = c("BMMC" = "gray", "PBMC" = "red")
)
SCpubr::do_ChordDiagramPlot(
 sample = nk,
 from = "condition",
 to = "seurat_clusters",
 colors.from = c("HD" = "green", "MM" = "purple")
)
SCpubr::do_ChordDiagramPlot(sample = nk, from = "orig.ident", to = "seurat_clusters")
SCpubr::do_ChordDiagramPlot(sample = nk, from = "id", to = "seurat_clusters")
