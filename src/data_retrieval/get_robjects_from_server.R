library(qs)
library(Seurat)


setwd(here::here())

###########################
###### Seurat object ######
###########################
seurat <- qread("/Volumes/CRCT13//Ruchan/01_NK_project/01_human/data/nonprolif_NKs_multiple_res_clusters_20dec22.qs")


seurat <- DietSeurat(seurat, assays = c("RNA", "ADT"), dimreducs = c("umap", "harmony"))

path_to_seurat_objects <- file.path("results","NK_analysis","seurat_objects")
dir.create(path_to_seurat_objects, recursive = T)
qsave(seurat, file.path(path_to_seurat_objects,"analysisNK_Ruchan.qs"))


###########################
#### slingshot object #####
###########################


# get slingshot results from Ruchan on the share
RuNK <- qread("/Volumes/CRCT13/Ruchan/01_NK_project/01_human/8_clusters_nonproliferating_nk/data/sling-test/density_st_c5_dont_allow_tree_break_sling.qs")

qsave(RuNK, "results/NK_analysis/trajectory_pseudotime/slingshot_pto.qs")

###########################
###### GSEA results #######
###########################



gsea_results <- qread(
  file.path(
    '/Volumes',
    'CRCT13',
    'Ruchan',
    '01_NK_project',
    '01_human',
    '8_clusters_nonproliferating_nk',
    'Enrichment_results',
    'detection_1_thresholded_auc_all_combin_gsea_param1_res.qs'
  )
)

qsave(gsea_results, "results/NK_analysis/clustering_typing/gsea_results.qs")