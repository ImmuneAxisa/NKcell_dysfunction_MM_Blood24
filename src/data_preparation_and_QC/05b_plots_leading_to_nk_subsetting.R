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
        library(ggplot2)
        library(ggalluvial)
        library(ggrepel)
        library(patchwork)
})
################################################################################
##    Starting with the loading of the end product of "05a_preprocessing.R"   ##
################################################################################

preprocess_plots <- "../pub_ready/docs/preprocess_plots/"
trim_sc <- qread(file = "../pub_ready/data/preprocessed_t_nk.qs")
meta_df <- trim_sc@meta.data %>% rownames_to_column(var = "cell_id")
# ==============================================

#' The granularity of the clusters are starting from resolution **0.1** and
#' reaching maximum as **1.0** with Leiden algorithm application,
#' our dataset consist on **22 unsupervised clusters**.
#'

#### Alluvial plot of cells' flow in different clustering resolutions ####
p <-
        SCpubr::do_AlluvialPlot(
                sample = trim_sc,
                first_group = "RNA_snn_res.0.1",
                middle_groups = c(
                        "RNA_snn_res.0.3",
                        "RNA_snn_res.0.5",
                        "RNA_snn_res.0.7",
                        "RNA_snn_res.0.9",
                        "RNA_snn_res.1"
                ),
                last_group = "predicted.celltype.l1",
                use_labels = T,
                font.size = 10,
                stratum.color = "black",
                repel = TRUE,
                stratum.fill.conditional = TRUE,
                legend.position = "bottom",
                use_viridis = FALSE,
                use_geom_flow = TRUE,
                curve_type = "sigmoid",
                plot.title = "Flow of the cells for different clusteing resolution, annotated by predicted cell types",
                plot.subtitle = "Numbers on strata are clusters, axis is clustering resolution except last axis of the predicted cell type level",
                colors.use = scCustomize::Hue_Pal(8) #"#F8766D" "#CD9600" "#7CAE00" "#00BE67" "#00BFC4" "#00A9FF" "#C77CFF" "#FF61CC"
        )

ggsave(
        filename = "Alluvial_t_nk_labelled_w_pred_cell_type1.pdf",
        p,
        device = "pdf",
        dpi = "retina",
        width = 16,height = 9,
        path = preprocess_plots
)

##### Clustering tree and labeling the most abundant "predicted cell type" of the given cluster ####

Mode <- function(x) {
        # ------ This function is extracting statistical mode of given data----------------
        # R-base mode() function is not useful, it gives mode of storage for given data.
        # This function is used within clustreee's "node_label_aggr" parameter to
        # describe how to assign aggregated cluster "label".
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
}
##
ct1 <- clustree::clustree(
                trim_sc,
                prefix = "RNA_snn_res.",
                assay = "RNA",
                exprs = "data",
                node_label = "predicted.celltype.l1",
                node_label_aggr = "Mode", # Node label parameter is requiring external function to call for "how to aggregate strings" per cluster
                node_label_size = 3
        ) + labs(subtitle = "Unsupervised clustering for different resolutions and 'mode' of predicted cell type labels")
ggsave(
        filename = "Clustree_t_nk_labelled_w_pred_cell_type1.pdf",
        ct1,
        device = "pdf",
        dpi = "retina",
        path = preprocess_plots,
        height = 10,
        width = 12
)
##
ct2 <-
        clustree::clustree(
                trim_sc,
                prefix = "RNA_snn_res.",
                assay = "RNA",
                exprs = "data",
                node_label = "predicted.celltype.l2",
                node_label_aggr = "Mode",# Node label parameter is requiring external function to call for "how to aggregate strings" per cluster
                node_label_size = 2
        ) + labs(subtitle = "Unsupervised clustering for different resolutions and 'mode' of predicted cell type labels 2")
ggsave(
        filename = "Clustree_t_nk_labelled_w_pred_cell_type2.pdf",
        ct2,
        device = "pdf",
        dpi = "retina",
        path = preprocess_plots,
        height = 10,
        width = 15
)

# =================================== #
##      Cluster  visualizations      ##
# =================================== #

# Representation of each cluster within donors
# annotated by predicted cell types level 1
p1 <- meta_df %>%
        ggplot(aes(id, fill = predicted.celltype.l1)) +
        geom_bar() +
        facet_wrap(~ seurat_clusters, scales = "free", ncol = 6) +
        theme(axis.text.x = element_text(
                angle = 90,
                vjust = 0.5,
                hjust = 0.5
        )) +
        labs(x = "Donor ID", y = "Cell Counts", title = "Per Cluster representation of the each donor")
ggsave(
        filename = "Cluster_representation_per_donor,pred_cell_type1.pdf",
        p1,
        device = "pdf",
        dpi = "retina",
        path = preprocess_plots,
        height = 10,
        width = 16
)
# Representation of each cluster within donors
# annotated by predicted cell types level 2
p2 <- meta_df %>%
        ggplot(aes(id, fill = predicted.celltype.l2)) +
        geom_bar() +
        facet_wrap(~ seurat_clusters, scales = "free", ncol = 6) +
        theme(axis.text.x = element_text(
                angle = 90,
                vjust = 0.5,
                hjust = 0.5
        )) +
        labs(x = "Donor ID", y = "Cell Counts", title = "Per Cluster representation of the each donor")
ggsave(
        filename = "Cluster_representation_per_donor,pred_cell_type2.pdf",
        p2,
        device = "pdf",
        dpi = "retina",
        path = preprocess_plots,
        height = 10,
        width = 18
)

p3 <- meta_df %>%
        ggplot(aes(x =id, y = percent.mt)) +
        geom_boxplot(outlier.size = 0) +
        facet_wrap(~ seurat_clusters, scales = "free", ncol = 6) +
        theme(axis.text.x = element_text(
                angle = 90,
                vjust = 0.5,
                hjust = 0.5
        )) +
        labs(x = "Donor ID", y = "Percentage of mitochondrial transcripts",
             title = "Per donor mitochondiral reads in each cluster")
ggsave(
        filename = "Per_donor_mito_reads_in_clusters.pdf",
        p3,
        device = "pdf",
        dpi = "retina",
        path = preprocess_plots,
        height = 10,
        width = 18
)
p4 <- meta_df %>%
        ggplot(aes(x =id, y = percent.rb)) +
        geom_boxplot(outlier.size = 0) +
        facet_wrap(~ seurat_clusters, scales = "free", ncol = 6) +
        theme(axis.text.x = element_text(
                angle = 90,
                vjust = 0.5,
                hjust = 0.5
        )) +
        labs(x = "Donor ID", y = "Percentage of ribosomal transcripts",
             title = "Per donor ribosomal reads in each cluster")
ggsave(
        filename = "Per_donor_ribo_reads_in_clusters.pdf",
        p4,
        device = "pdf",
        dpi = "retina",
        path = preprocess_plots,
        height = 10,
        width = 18
)

#'
#' For the sake of better distinction of NK cells from the T cells, unsupervised clusters coming from the granularity of 1.0 is used.
#'

u1 <-
        SCpubr::do_DimPlot(
                trim_sc,
                reduction = "umap",
                group.by = "seurat_clusters",
                shuffle = TRUE,
                pt.size = 1,
                legend.position = "none",
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
ggsave(
        filename = "Unsupervised_clustering_t_nk_comb_umap.pdf",
        u1,
        device = "pdf",
        dpi = "retina",
        path = preprocess_plots,
        height = 10,
        width = 18
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
                trim_sc,
                reduction = "umap",
                group.by = "predicted.celltype.l1",
                shuffle = TRUE,
                pt.size = 1,
                legend.title = "Inferred cell type annotations",legend.position = "right",
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
ggsave(
        filename = "Predicted_main_cell_types_t_nk_comb_umap.pdf",
        u2,
        device = "pdf",
        dpi = "retina",
        path = preprocess_plots,
        height = 10,
        width = 18
)

u3 <-
        SCpubr::do_DimPlot(
                trim_sc,
                reduction = "umap",
                group.by = "predicted.celltype.l2",
                shuffle = TRUE,
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
ggsave(
        filename = "Predicted_cell_subtypes_t_nk_comb_umap.pdf",
        u3,
        device = "pdf",
        dpi = "retina",
        path = preprocess_plots,
        height = 10,
        width = 18
)

#' On the other hand, there was one donor (T23491) which undergone sequencing of T cells and NK cells, as separate emulsions.\
#' Usage of those distinct runs can be traced and we can use those cells to assess our unsupervised investigation of NK cells, like a positive control.\
#' This donor is practically a ground truth for our cell type classification.\
#' Below are the figures of that donor with different "perspectives".

#''------------------------------------------------------------------------------------------------------------------
c1 <- SCpubr::do_DimPlot(sample = trim_sc[,trim_sc$orig.ident == "T23491b_NK"], #This sample is expected to contain only NK cells
        reduction = "umap",
        group.by = "predicted.celltype.l2",
        shuffle = TRUE,
        pt.size = 1,
        legend.title = "Inferred Annotation",legend.position = "none",
        legend.icon.size = 3,font.size = 12,font.type = "serif",
        plot.axes = FALSE,plot.title = "T23491b_NK",
        plot_density_contour = T,
        contour.color = "gray",
        plot_cell_borders = T,
        raster = T,
        contour.position = "top",
        label = T,
        repel = TRUE,
        label.box = TRUE
)
c2 <- SCpubr::do_DimPlot(sample = trim_sc[,trim_sc$orig.ident == "T23491b_T"], #This sample is expected to contain only T cells
                         reduction = "umap",
                         group.by = "predicted.celltype.l2",
                         shuffle = TRUE,
                         pt.size = 1,
                         legend.title = "Inferred Annotation",legend.position = "none",
                         legend.icon.size = 3,font.size = 12,font.type = "serif",
                         plot.axes = FALSE,plot.title = "T23491b_T",
                         plot_density_contour = T,
                         contour.color = "gray",
                         plot_cell_borders = T,
                         raster = T,
                         contour.position = "top",
                         label = T,
                         repel = TRUE,
                         label.box = TRUE
)
c3 <- SCpubr::do_DimPlot(sample = trim_sc[,trim_sc$orig.ident == "T23491m_NK"], #This sample is expected to contain only NK cells
                         reduction = "umap",
                         group.by = "predicted.celltype.l2",
                         shuffle = TRUE,
                         pt.size = 1,
                         legend.title = "Inferred Annotation",legend.position = "none",
                         legend.icon.size = 3,font.size = 12,font.type = "serif",
                         plot.axes = FALSE,plot.title = "T23491m_NK",
                         plot_density_contour = T,
                         contour.color = "gray",
                         plot_cell_borders = T,
                         raster = T,
                         contour.position = "top",
                         label = T,
                         repel = TRUE,
                         label.box = TRUE
)
c4 <- SCpubr::do_DimPlot(sample = trim_sc[,trim_sc$orig.ident == "T23491m_T"], #This sample is expected to contain only T cells
                         reduction = "umap",
                         group.by = "predicted.celltype.l2",
                         shuffle = TRUE,
                         pt.size = 1,
                         legend.title = "Inferred Annotation",legend.position = "none",
                         legend.icon.size = 3,font.size = 12,font.type = "serif",
                         plot.axes = FALSE,plot.title = "T23491m_T",
                         plot_density_contour = T,
                         contour.color = "gray",
                         plot_cell_borders = T,
                         raster = T,
                         contour.position = "top",
                         label = T,
                         repel = TRUE,
                         label.box = TRUE
)
# Patchwork of these plots
(c1+c2)/(c3+c4)

# Barplot of the positive control samples' inferred cell types
b1 <-   meta_df %>% dplyr::filter(id == "T23491") %>%
        complete(orig.ident, predicted.celltype.l1, fill = list(predicted.celltype.l2 = NA, count = 0)) %>%
        ggplot(aes(y = predicted.celltype.l1, fill = predicted.celltype.l2)) +
        geom_bar(position = "stack") + scale_fill_discrete(na.translate=FALSE) +
        facet_wrap(
                ~ orig.ident,
                scales = "free",drop = FALSE,
                ncol = 2
        ) +
        theme(axis.text.x = element_text(
                angle = 90,
                vjust = 0.5,
                hjust = 0.5
        )) +
        labs(
                x = "Cell counts",
                y = "Cell types",
                # title = "Per cell type representation of the 'positive control' donor",
                fill = "Inferred cell sub-types"
        )
ggsave(
        filename = "Cell_types_of_pos_control_barplot.pdf",
        b1,
        device = "pdf",
        dpi = "retina",
        path = preprocess_plots,
        height = 6,
        width = 9
)
#' Lets look at marker genes of these clusters.
#'
## Cluster markers
cluster_markers <-
        FindAllMarkers(
                trim_sc,
                assay = "RNA",
                test.use = "MAST",
                min.pct = 0.25,
                logfc.threshold = 0.25,
                verbose = T
        )
qsave(cluster_markers,file = "data/clusterMarkersTPlusNK.qs")

# Access the curated gene symbol annotations
annotations <- qs::qread(file = "data/gene_annotations_of_this_run.qs")
#'
cluster_markers %>% scCustomize::Add_Pct_Diff(overwrite = T) %>%
        group_by(cluster) %>%
        top_n(n = 10, wt = avg_log2FC) -> top10_cluster_marker
top10_cluster_marker %>%
        left_join(
                x = top10_cluster_marker,
                y = unique(annotations[, c("gene_name", "description")]),
                by = c("gene" = "gene_name")
        )
top10_cluster_marker


################
genes_to_plot <- trim_sc@assays$RNA@data[unique(top5_cluster_marker$gene),] %>% #marker_genes_to_plot
 as.matrix() %>%
 t() %>%
 as.data.frame() %>%
 dplyr::rename_with(~paste0("gene_",.)) %>%
 tibble::rownames_to_column("cell.id")

meta <- trim_sc@meta.data

hm_genes <- meta %>%
 tibble::rownames_to_column("cell.id") %>%
 dplyr::select(cell.id, seurat_clusters) %>%
 left_join(genes_to_plot) %>%
 gather("feat", "val",  -!matches("^gene_")) %>%
 group_by(feat, seurat_clusters) %>%
 summarize(avg = mean(val))


hm_genes <- hm_genes %>%
 mutate(feat = gsub("gene_", "", feat)) %>%
 pivot_wider(names_from = "feat",values_from = "avg") %>%
 tibble::column_to_rownames("seurat_clusters") %>%
 as.matrix() %>%
 scale()

pdf("averageExptTop5MarkerHeatmapTPLusNKlog2fcMast.pdf",paper = "a4")
library(ComplexHeatmap)
hm <- ComplexHeatmap::Heatmap(t(hm_genes),
                              col = circlize::colorRamp2(quantile(hm_genes, c(0.025,0.5,0.975)),
                                               c("dodgerblue3", "white", "firebrick3")),name = "Relative gene expression",
                              row_names_gp = grid::gpar(fontsize = 4, fontface = "bold"),height = 16,
                              column_names_gp = grid::gpar(fontsize = 8, fontface = "bold"),
                              column_names_rot = 0,cluster_columns = T,clustering_distance_columns = "pearson",
                              cluster_rows = T,clustering_distance_rows = "pearson",
                              column_names_centered = T,row_names_centered = T,use_raster = F,
                              column_split = 4,row_split = 4,heatmap_legend_param = list(
                               legend_direction = "horizontal",just = "center",
                               legend_width = unit(10, "cm")))

hm <- draw(hm,heatmap_legend_side="bottom")
dev.off()

#' ## UMAPs of Interesting Features
#' ### Protein expressions
## ----results='hide'----------------------------------------------------------------------------------------------------------------------------------------------
DefaultAssay(trim_sc) <- "ADT"
# trim_sc$ident<- NULL
Idents(trim_sc) <- trim_sc$seurat_clusters
FeaturePlot(
        trim_sc,
        features = sort(c(
                rownames(trim_sc@assays$ADT)[1:9],
                rownames(trim_sc@assays$ADT)[12:14]
        )),
        order = T,
        label = T,
        ncol = 4,
        combine = F,
        cols = scCustomize::viridis_magma_dark_high
)
DefaultAssay(trim_sc) <- "RNA"


#### Protein exprs to plot ####
prots_to_plot <- trim_sc@assays$ADT@data %>% #marker_genes_to_plot
 as.matrix() %>%
 t() %>%
 as.data.frame() %>%
 dplyr::rename_with(~paste0("gene_",.)) %>%
 tibble::rownames_to_column("cell.id")

meta <- trim_sc@meta.data

hm_genes <- meta %>%
 tibble::rownames_to_column("cell.id") %>%
 dplyr::select(cell.id, seurat_clusters) %>%
 left_join(prots_to_plot) %>%
 gather("feat", "val",  -!matches("^gene_")) %>%
 group_by(feat, seurat_clusters) %>%
 summarize(avg = mean(val))


hm_genes <- hm_genes %>%
 mutate(feat = gsub("gene_", "", feat)) %>%
 pivot_wider(names_from = "feat",values_from = "avg") %>%
 tibble::column_to_rownames("seurat_clusters") %>%
 as.matrix() %>%
 scale()

pdf("averageExpProtHeatmapTPLusNK.pdf",paper = "a4")
library(ComplexHeatmap)
hm <- ComplexHeatmap::Heatmap(t(hm_genes),
                              col = circlize::colorRamp2(quantile(hm_genes, c(0.025,0.5,0.975)),
                                                         c("dodgerblue3", "white", "firebrick3")),name = "Relative cell-surface protein expression",
                              row_names_gp = grid::gpar(fontsize = 10, fontface = "bold"),height = 16,
                              column_names_gp = grid::gpar(fontsize =8 , fontface = "bold"),
                              column_names_rot = 0,cluster_columns = T,clustering_distance_columns = "pearson",
                              cluster_rows = T,clustering_distance_rows = "pearson",
                              column_names_centered = T,row_names_centered = T,use_raster = F,
                              column_split = 4,row_split = 4,heatmap_legend_param = list(
                               legend_direction = "horizontal",just = "center",
                               legend_width = unit(8, "cm")))

hm <- draw(hm,heatmap_legend_side="bottom")
dev.off()



#' ### Gene expressions
## ----results='hide'----------------------------------------------------------------------------------------------------------------------------------------------
trim_sc$ident <- trim_sc$seurat_clusters
FeaturePlot(
        trim_sc,
        features = sort(
                c(
                        "GZMB",
                        "ZBTB16",
                        "GNLY",
                        "NKG7",
                        "KLRD1",
                        "FCER1G",
                        "SPON2",
                        "CD4",
                        "FCGR3A",
                        "PRF1",
                        "CD8A",
                        "CD8B",
                        "CD3E",
                        "CD3D",
                        "CD3G"
                )
        ),
        order = T,
        label = T,
        ncol = 4,
        combine = F,
        cols = scCustomize::viridis_magma_dark_high
)


#'
#'
#'
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
