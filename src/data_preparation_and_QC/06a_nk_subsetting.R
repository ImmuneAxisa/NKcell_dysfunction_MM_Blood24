#' title: "Subsetting NK cells"
##########################
##   Loading packages   ##
##########################
suppressPackageStartupMessages({
  library(tidyverse)
  library(conflicted)
  library(qs)
  library(Seurat)
  library(tidyseurat)
  library(harmony)
  library(clustree)
  library(ggplot2)
  library(ggalluvial)
  library(ggrepel)
  library(patchwork)
})
#############################################################################
##  After seeing the convincing evidences of "predicted cell-type assignment
##  is pretty accurate for NK cells, we decided to subset all NK main-type
##  assigned cells for this part.-------------------------------------------
##  Starting with the loading of the end product of "05a_preprocess.R"  ##
############################################################################
# Making it ready for the subsetting of the NK cells.
preprocess_plots <- "../pub_ready/docs/preprocess_plots/"
trim_sc <- qread(file = "../pub_ready/data/preprocessed_t_nk.qs")

## ----NK subsetting------------------------------------------------------------
nk_part <- subset(trim_sc,predicted.celltype.l1 == "NK")

# ---- Discarding residual bits, carried from before subsetted dataset----------
nk_part <-DietSeurat(object = nk_part,counts = T,scale.data = F,dimreducs = NULL,graphs = NULL)
nk_part

nk_part %>% tidyseurat::as_tibble() %>% select(!contains("_res")) %>% select(!starts_with("scDblFinder"))
#'
#' Choosing justifiable \# dims of Harmony reduction,
#' It is highlighted that first 29 harmony planes
#' (dimensions) are optimal to be worked with.
## --------------------------------------------------------------
nk_part <-
  nk_part %>% NormalizeData() %>% FindVariableFeatures() %>%
  ScaleData(vars.to.regress = c("percent.mt", "percent.rb"))  %>%
  RunPCA()
nk_part <- nk_part %>% RunHarmony(group.by.vars = "id")

source(file = "../pub_ready/func/quant_pcs.R");source(file = "../pub_ready/func/quant_harmonys.R")
quant_opt_pcs(nk_part); quant_opt_harmonys(nk_part)

#'---------------------------------------------------------------------------------------------------------------------------
# First 22 harmony corrected PC dimensions are chosen for neighborhood detection and clustering
#'--------------------------------------------------------------------------------------------------------------------------
nk_part <- FindNeighbors(nk_part,reduction = "harmony",dims = 1:22,force.recalc = T)
nk_part <- FindClusters(object = nk_part,algorithm = 4,method = "igraph",resolution = seq(0.1,1,0.1),random.seed = 2206,verbose = T)
nk_part <- RunUMAP(nk_part, reduction = "harmony", dims = 1:22, umap.method = "umap-learn", metric = "correlation")
qsave(nk_part,file = "allNKsubsetsMultiResolutionClustered.qs")
## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
u1 <-
  SCpubr::do_DimPlot(
    nk_part,
    reduction = "umap",
    group.by = "seurat_clusters",
    shuffle = F,
    pt.size = 1,
    legend.position = "none",
    plot.axes = FALSE,
    plot_density_contour = T,
    contour.color = "gray",
    plot_cell_borders = T,
    raster =T,
    contour.position = "top",
    label = T,
    repel = TRUE,
    label.box = TRUE
  )
ggsave(
  filename = "Unsupervised_clustering_nk_part_umap.pdf",
  u1,
  device = "pdf",
  dpi = "retina",
  path = preprocess_plots,
  height = 8,
  width = 14
)
#
# Apart from unsupervised clustering approach,
# I have "label transfer"ed the whole PBMC CITEseq reference data to our dataset in the very beginning of the analysis
#
#' Below figures are the visualization of predicted cell types within our dataset.
#'
## ----fig.height=8, fig.width=10, results='hide'------------------------------------------------------------------------------------------------------------------
u2 <-
  SCpubr::do_DimPlot(
    nk_part,
    reduction = "umap",
    group.by = "predicted.celltype.l2",
    shuffle = F,
    pt.size = 1,
    legend.title = "Inferred cell sub-type annotations",legend.position = "right",
    legend.icon.size = 3,font.size = 12,font.type = "serif",
    plot.axes = FALSE,
    plot_density_contour = T,
    contour.color = "gray",
    plot_cell_borders = T,
    raster = T,
    contour.position = "top",
    label = T,
    repel = TRUE,
    label.box = TRUE
  )
u1+u2
# qsave(nk_part,file = "nk_part_clustered")
nk_part <- qread(file = "data/denoised_clustered_nk",nthreads = 10)

# ------------------------------------------------------------------------------------------------------------------------------
all_nk_markers <- FindAllMarkers(nk_part,assay = "RNA",test.use = "MAST",min.pct = 0.25, logfc.threshold = 0.25,verbose = T,only.pos = T)
#'
## ----gene annotations-------------------------------------------------------------------------------------------------------------------------------------

# Select annotations of interest
annotations <- annotations %>%
        dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)
# -------------------------------------------------------------------------------------------------------------------------------------------------------
nk_part <- RunUMAP(nk_part,reduction = "harmony",dims = 1:22,umap.method = "umap-learn",metric = "correlation",min.dist = 0.01,spread = 1)
# qsave(nk_part,file = "denoised_clustered_nk")
nk_part <- qread(file = "denoised_clustered_nk")

#'
#' ### Visualizations
#'
#'
#'------------------------------------------------------------------------------------------------------------------------------


#'----------------------------------------------------------------------------------------------------------------------------
VlnPlot(object = nk_part,
        features = sort(unique(c("NKG7","GNLY","FCGR3A","NCAM1","PRF1","KLRB1","KLRD1","KLRF1","KLRC2","CX3CR1", "CXCR4", "CXCR6","GZMA","GZMB","GZMK","TIGIT","LAG3","FCER1G", "CXCR4","KLRC1","MAPK3","ITGAL","SELL","CCR7","SPON2","CD2","MAPK3","ITGAL","ZBTB16","CD3E","CD3D","IL2","TNF","IFNG","IL4","CD2","CD3G","CD4","CD8A","CD8A","CD8B","CD247"))),stack = T,flip = T) + NoLegend()

#'
#' ### Markers of clusters
# ------------------------------------------------------------------------------
all_nk_markers  %>% scCustomize::Add_Pct_Diff(overwrite = T) %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC) -> top30_nk_cluster_marker

#'
#' Top 30 DEGs by the avg log2FC of the genes compared to other clusters.
#' *pct.1* indicates how much percent of the cells on this cluster positively express the gene of interest, and *pct.2* is the averaged percent expression that gene in other clusters.
# ------------------------------------------------------------------------------
top30_nk_cluster_marker%>% left_join(
  x = top30_nk_cluster_marker[,-1],
  y = unique(annotations[, c("gene_name", "description")]),
  by = c("gene" = "gene_name")
)

#' Top 30 markers by weighted presence on a cluster
# ------------------------------------------------------------------------------
all_nk_markers  %>% scCustomize::Add_Pct_Diff(overwrite = T) %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = pct_diff) -> top30_exc_cluster_marker
top30_exc_cluster_marker %>% left_join(
  x = top30_exc_cluster_marker[,-1],
  y = unique(annotations[, c("gene_name", "description")]),
  by = c("gene" = "gene_name")
)
## ----Proliferating NK Genes, eval=FALSE, include=FALSE--------------------------------------------------------------------------------------------------------------------
## VlnPlot(nk_part,features = c("STMN1","MKI67","CENPF","PCLAF","TYMS","HMGN2","PCNA"),pt.size = 0,group.by = "seurat_clusters",slot = "data",stack = T,flip = T)+NoLegend()
#'
#' ### UMAPs of Interesting Features
#' #### Protein expressions
#'-------------------------------------------------------------------------------------------------------------------------------------------------------
DefaultAssay(nk_part)<-"ADT"
nk_part$ident<- NULL
Idents(nk_part) <- nk_part$seurat_clusters
FeaturePlot(nk_part,features = sort(c(rownames(nk_part@assays$ADT)[1:9],rownames(nk_part@assays$ADT)[12:14])),
            label=T,ncol = 4,combine = F,cols = scCustomize::viridis_magma_dark_high)
DefaultAssay(nk_part)<-"RNA"

#'
#' #### Gene expressions
##-------------------------------------------------------------------------------------------------------------------------------------------------------

nk_part$ident<- nk_part$seurat_clusters
FeaturePlot(nk_part,features = sort(unique(c("NKG7","GNLY","FCGR3A","NCAM1","PRF1","KLRB1","KLRD1","KLRF1","KLRC2","CX3CR1", "CXCR4", "CXCR6","GZMA","GZMB","GZMK","TIGIT","LAG3","FCER1G", "CXCR4","KLRC1","MAPK3","ITGAL","SELL","CCR7","SPON2", "CD2","MAPK3","ITGAL","ZBTB16","CD3E","CD3D","IL2","TNF","IFNG","IL4","CD2","CD3G","CD4","CD8A","CD8A","CD8B","CD247"))),
            label = T,ncol = 4,combine = F,cols = scCustomize::viridis_magma_dark_high)


