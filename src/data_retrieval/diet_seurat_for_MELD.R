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
in_stem <- "analysisNK_Ruchan"

input_seurat <- file.path(objects_dir,paste0(in_stem,".qs"))

#--------------------
# output path

#output name stem
out_stem <- in_stem

out_seurat <- file.path(
  "results",
  "NK_analysis",
  "seurat_objects",
  paste0("seurat_diet_",out_stem,".rds")
)

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



#########################################################
###################### import data ######################
#########################################################


# get seurat object
merged.seurat <- qread(input_seurat)

meta <- merged.seurat@meta.data %>%
  tibble::rownames_to_column("cell.id")

merged.seurat@meta.data <- meta %>%
  curation() %>%
  tibble::column_to_rownames("cell.id")


merged.seurat <- DietSeurat(merged.seurat, dimreducs = "harmony")


#########################################################
###################### export data ######################
#########################################################

saveRDS(merged.seurat,
        out_seurat)

