#################################################
################### Set up ######################
#################################################

## load libraries
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(data.table)

# export path
export_csv_dir <- "results/NK_analysis/phate_meld_input_csv/"

#########################################################
############ write pca and metadata table  ##############
#########################################################

#-----------
# import data
#-----------
merged.seurat <- readRDS("results/NK_analysis/seurat_objects/seurat_diet_analysisNK_Ruchan.rds")


#---------------------------
# merge pc loadings and meta
#---------------------------
master_harmony <- merged.seurat@reductions$harmony@cell.embeddings %>% 
		as.data.frame() %>%
		tibble::rownames_to_column("cell.id") %>%
		left_join(merged.seurat@meta.data %>% tibble::rownames_to_column("cell.id"), by = "cell.id") %>%
        dplyr::rename("index" = "cell.id") %>%
        dplyr::rename("tissue" = "source")
        


#-----------------------------        
# define splitting categories
#-----------------------------

tissue_cats <- unique(master_harmony$tissue)
condition_cats <- unique(master_harmony$condition)


#--------------------------------------------  
# filter on tissue and run MELD on disease
#--------------------------------------------  
for(tissue_index in tissue_cats) {


    # do filterting and cleaning before export for phate meld vfc        
            
    harmony <- master_harmony %>%
            dplyr::filter(tissue %in% tissue_index)
    
    
    # split and write to file
      
    harmony_PCA <- harmony[,1:31]
    
    harmony_meta <- harmony[,c(1,52:ncol(harmony))]
    
    basename <- paste(condition_cats, collapse = "vs")
    basename <- paste(basename, tissue_index, sep = "_")
    
    write_csv(harmony_PCA, paste0(export_csv_dir,basename,"_pca.csv"))
    
    write_csv(harmony_meta, paste0(export_csv_dir,basename,"_meta.csv"))
    
    
    # run the shell wrapper to submit phate and meld jobs
    
    
    system2(command = "src/MELD/phate_meld_sbatcher_w_args_v2.sh", 
            args    = c("condition", #column to use in meta to run MELD with
                        "curated_cluster", #column to use in meta to split data into for the VFC 
                        "NK", #suffix to add for the output
                        basename
                        ))
}


#--------------------------------------------
# filter on disease and run MELD on tissue
#--------------------------------------------
for(condition_index in condition_cats) {


    # do filterting and cleaning before export for phate meld vfc        
            
    harmony <- master_harmony %>%
            dplyr::filter(condition %in% condition_index)
    
    
    # split and write to file
      
    harmony_PCA <- harmony[,1:31]
    
    harmony_meta <- harmony[,c(1,52:ncol(harmony))]
    
    basename <- paste(tissue_cats, collapse = "vs")
    basename <- paste(basename, condition_index, sep = "_")
    
    write_csv(harmony_PCA, paste0(export_csv_dir,basename,"_pca.csv"))
    
    write_csv(harmony_meta, paste0(export_csv_dir,basename,"_meta.csv"))
    
    
    # run the shell wrapper to submit phate and meld jobs
    
    
    system2(command = "src/MELD/phate_meld_sbatcher_w_args_v2.sh", 
            args    = c("tissue", #column to use in meta to run MELD with
                        "curated_cluster", #column to use in meta to split data into for the VFC 
                        "NK", #suffix to add for the output
                        basename
                        ))
}

