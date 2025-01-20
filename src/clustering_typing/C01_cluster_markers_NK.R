#########################################################
################# Input/output info #####################
#########################################################

setwd(here::here())

#--------------------
# directory structure
results_dir <- "results/NK_analysis"
dir.create(results_dir, showWarnings = FALSE)

objects_dir <- file.path(results_dir,"seurat_objects")
dir.create(objects_dir, showWarnings = FALSE)

clust_dir <- file.path(results_dir,"clustering_typing")
dir.create(clust_dir, showWarnings = FALSE)


#--------------------
# input path

input_seurat <- file.path(objects_dir, "analysisNK_Ruchan.qs")

#input_table <- file.path(clust_dir,"clustering_table_harmony_CD8T_regressRb.qs")

#--------------------
# output path

#output name stem
out_stem <- "NK"

# table output
output_table <- file.path(clust_dir,paste0("clusters_markers_table_harmony_", out_stem,".qs"))

output_table_filtered <- file.path(clust_dir,paste0("clusters_markers_table_harmony_", out_stem,"_filtered.qs"))

output_table_paired <- file.path(clust_dir,paste0("clusters_markers_table_harmony_", out_stem,"_paired.qs"))

output_table_paired_filtered <- file.path(clust_dir,paste0("clusters_markers_table_harmony_", out_stem,"_paired_filtered.qs"))


#--------------------
# options


# function to curate cluster calls to test: needs to generate a "curated_cluster" column in the meta
curation <- function(meta) {
  meta %>%
    mutate(curated_cluster = seurat_clusters) 
  
}

#################################################
################### Set up ######################
#################################################

## Take the input argument
# args <- commandArgs(TRUE)
# sample_name <- args[1]
# exp_matrix <- args[2]
# demux_best <- args[3]
# out_path <- args[4]

## load libraries

library(gridExtra)
library(knitr)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggExtra)
library(stats)
# library(RColorBrewer)
library(viridis)
library(scales)
library(ggrepel)
library(circlize)
library(RColorBrewer)
library(gplots)
library(cowplot)
library(gtable)
library(data.table)
library(Seurat)
library(hdf5r)
library(ggsci)
library(harmony)
library(presto)
library(qs)



## ggplot theme

theme_mrl <- function(x = 1) {
  theme_minimal() +
    theme(
      axis.line = element_line(),
      axis.ticks.x = element_line(),
      axis.ticks.y = element_line(),
      axis.text.x = element_text(
        size = 12 * x,
        face = "bold",
        angle = 45,
        vjust = 0.6
      ),
      axis.text.y = element_text(size = 12 * x, face = "bold"),
      axis.title.x = element_text(size = 12 * x, face = "bold"),
      axis.title.y = element_text(size = 12 * x, face = "bold"),
      strip.background = element_rect(
        fill = "gray20",
        colour = "gray20",
        linetype = "solid"
      ),
      strip.text = element_text(
        size = 14 * x,
        colour = "white",
        face = "bold"
      ),
      legend.title = element_text(size = 14 * x, face = "bold"),
      legend.text = element_text(
        size = 12 * x,
        color = "gray20",
        face = "bold"
      ),
      legend.background = element_rect(fill = "transparent", colour = "transparent"),
      plot.title =  element_text(
        hjust = 0.5,
        vjust = 2,
        face = "bold"
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        vjust = 3,
        face = "italic"
      ),
      plot.caption = element_text(hjust = 0, face = "italic")
    )
}





#########################################################
###################### import data ######################
#########################################################


# get seurat object
merged.seurat <- qread(input_seurat)


meta <- merged.seurat@meta.data %>%
  tibble::rownames_to_column("cell.id") %>%
  curation() %>%
  tibble::column_to_rownames("cell.id")

merged.seurat@meta.data <- meta

############################
####### find markers #######
############################


#----------------------------
# general DGE for each cluster


presto_results <-
  wilcoxauc(merged.seurat, "curated_cluster", assay = 'data')

qsave(
  presto_results,
  output_table
)

presto_results <- presto_results %>%
  dplyr::filter(padj < 0.01 & logFC > 0 & auc > 0.6)

qsave(
  presto_results,
  output_table_filtered
)


#----------------------------
# DGE paired for each cluster


cluster_assignments <- merged.seurat@meta.data$curated_cluster

combinations <- combn(unique(cluster_assignments),
                      2,
                      simplify = F)

paired_results <- lapply(combinations, function(comb) {
  wilcoxauc(merged.seurat,
            "curated_cluster",
            groups_use = comb,
            assay = 'data') %>%
    mutate(groupA = comb[1]) %>%
    mutate(groupB = comb[2])
})

paired_results <- do.call(rbind, paired_results)

qsave(
  paired_results,
  output_table_paired
)

paired_results <- paired_results %>%
  dplyr::filter(padj < 0.01 & logFC > 0 & auc > 0.6)

qsave(
  paired_results,
  output_table_paired_filtered
)



