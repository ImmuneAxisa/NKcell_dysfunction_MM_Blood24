###### Loading required packages ########
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(SeuratDisk)
  library(SeuratWrappers)
  library(BiocParallel)
  library(scDblFinder)
  library(qs)
})

#### File/Data Path Preps ####
# We downloaded the multimodal reference pbmc dataset from here  (https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat)
# This dataset is published in this article https://doi.org/10.1016/j.cell.2021.04.048
#
ref_pbmc_path <- "/mnt/SERVER-CRCT-STORAGE/CRCT13/Ruchan/01_NK_project/99_back_up/reference_databases/pbmc_multimodal_reference.h5seurat"
sample_paths <-  "/mnt/SERVER-CRCT-STORAGE/CRCT13/Ruchan/01_NK_project/99_back_up/cellranger_counts_human/"

# Loading up the reference pbmc dataset
ref_pbmc <- LoadH5Seurat(file = ref_pbmc_path)
Idents(ref_pbmc) <- ref_pbmc$celltype.l1

# Reading the samples of our cohort
name_of_samples <- list.files(sample_paths)
#### Main Task ####
# Handling each sample individually for ref-based cell-type annotation,
# doublet detection, velocity embedding and saving the processed file accodingly
for (i in 1:length(name_of_samples)) {
  ## Preprocessing steps
  filt.matrix <- Read10X(data.dir = paste0(sample_paths, name_of_samples[i],"/outs/filtered_feature_bc_matrix/"),strip.suffix = T)
  tmp <- CreateSeuratObject(counts = filt.matrix$`Gene Expression`, project = name_of_samples[i])
  tmp <- NormalizeData(tmp) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(reduction.name = "pca", verbose = F) %>%
    RunUMAP(dims = 1:20, verbose = F) %>% FindNeighbors(dims = 1:20, verbose = F) %>%
    FindClusters(algorithm = 3, verbose = T)
  ####
  #### ScDblFinder is a package generating artificial doublets by given clusters.
  #### and then training a model to classify singlets and doublets.
  #### https://f1000research.com/articles/10-979
  ####
  tmp_sce <- as.SingleCellExperiment(tmp) #ScdblFinder is compatible within SingleCellExperiment ecosysten!;
  set.seed(22)
  tmp_sce <-
    scDblFinder(tmp_sce,
                clusters = tmp_sce$seurat_clusters,
                BPPARAM = MulticoreParam(16))
  table(tmp_sce$scDblFinder.class)
  rm(tmp); gc()
  tmp <-
    as.Seurat(tmp_sce, project = (name_of_samples[i]))
  rm(tmp_sce); gc()

  #### Metadata ####
  tmp$id <- ifelse(tmp$orig.ident %in% c("T23491m_NK", "T23491m_T", "T23491b_NK", "T23491b_T"),
                   "T23491", str_sub(name_of_samples[i],1,nchar(name_of_samples[i])-1))
  tmp$percent.mt <- PercentageFeatureSet(tmp, pattern = "^MT-") # calculate and add mito read percentage
  tmp$percent.rb <- PercentageFeatureSet(tmp, pattern = "^RP[SL]") # calculate and add ribo read precentage
  # If the donor was Healthy, sample name was starting with "H",
  #  else sample name was starting with "T"
  # HD: Healthy Donor, MM: Multiple Myeloma
  tmp$condition = ifelse(grepl(pattern="H", x = name_of_samples[i]), "HD", "MM")
  # If the sample was taken from Bone marrow, it was encoded at the end of sample name with "m",
  # else name was ending with "b" which reflects peripheral blood
  # BMMC: Bone Marrow Mononuclear Cells, PBMC: Peripheral Blood Mononuclear Cells
  tmp$source = ifelse(grepl(pattern="m", x = name_of_samples[i]), "BMMC", "PBMC") #
  ####
  DefaultAssay(tmp) <- "RNA"
  tmp@reductions$PCA <- NULL
  tmp <-
    NormalizeData(tmp) %>% FindVariableFeatures() %>% SCTransform(
      method = "glmGamPoi",
      vst.flavor = "v2",
      vars.to.regress = "percent.mt",
      assay = "RNA",
      do.scale = T
    ) %>% RunPCA(reduction.name = "pca")
  tmp[["ADT"]] <-
    CreateAssayObject(filt.matrix$`Antibody Capture`[, colnames(x = tmp)])
  DefaultAssay(tmp) <- "ADT"
  VariableFeatures(tmp) <- rownames(tmp[["ADT"]])
  tmp <-
    NormalizeData(tmp, normalization.method = "CLR", margin = 2) %>%
    ScaleData() %>% RunPCA(reduction.name = "apca") ## Normalize protein data across cells
  DefaultAssay(tmp) <- "SCT"
  #Due to the characteristic of the used reference dataset,
  #cell type annotation accuracy is getting better with multimodal cell type label transfer approach.

  # To identify multimodal neighbors. These will be stored in the neighbors slot,
  # and can be accessed using sc_obj[['weighted.nn']]
  # The WNN graph can be accessed at sc_obj[["wknn"]],
  # and the SNN graph used for clustering at sc_obj[["wsnn"]]
  # Cell-specific modality weights can be accessed at sc_obj$RNA.weight
  tmp <- FindMultiModalNeighbors(tmp,
                                 reduction.list = list("pca", "apca"),
                                 dims.list = list(1:20, 1:10))

  tmp <-
    RunUMAP(
      tmp,
      nn.name = "weighted.nn",
      reduction.name = "wnn.umap",
      reduction.key = "wnnUMAP_"
    )
  tmp <-
    FindClusters(tmp,
                 graph.name = "wsnn",
                 algorithm = 3,
                 verbose = FALSE)


  anchors <- FindTransferAnchors(
    reference = ref_pbmc,
    query = tmp,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50
  )


  tmp <- TransferData(
    anchorset = anchors,
    reference = ref_pbmc,
    query = tmp,
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      predicted_ADT = "ADT"
    )
  )

  head(colnames(tmp))
  tmp <-
    RenameCells(object = tmp,
                new.names = paste0(name_of_samples[i], ":", colnames(x = tmp)))
  head(colnames(tmp))
  rm(anchors)
  gc()

  # Reading calculated RNA velocity information into the data object
  # Velocity assay embeddings, echo=TRUE, message=FALSE, warning=FALSE, paged.print=TRUE}
  tmp_velo <-
    ReadVelocity(file = paste0(sample_paths, name_of_samples[i], "/velocyto/", name_of_samples[i], ".loom"))
  tmp_velo <- as.Seurat(tmp_velo)
  tmp_velo$orig.ident <- as.factor(name_of_samples[i])
  #This "Renaming the cells" for uniform formatting is required to merge two object by their common cell ids.
  # Gene expression Seurat Object and Velocity Seurat objects can be merged thanks to this. #
  tmp_velo <- RenameCells(object = tmp_velo,
                          new.names = substring(colnames(tmp_velo), 1, nchar(colnames(tmp_velo)) -
                                                  1))
  tmp[["spliced"]] <-
    CreateAssayObject(tmp_velo[["spliced"]][, colnames(x = tmp_velo)])
  tmp[["unspliced"]] <-
    CreateAssayObject(tmp_velo[["unspliced"]][, colnames(x = tmp_velo)])
  tmp[["ambiguous"]] <-
    CreateAssayObject(tmp_velo[["ambiguous"]][, colnames(x = tmp_velo)])
  rm(tmp_velo)
  gc()
  print(name_of_samples[i])

  SaveH5Seurat(object = tmp, filename = paste0("/data/velo_embedded_wnn/",name_of_samples[i]))
  rm(tmp)
  gc()


}

#### Footnote ####
# Next step will be inference of Nuclear RNA fraction from Velocyto CLI outputs, directly.
# Then reducing this metric into a metadata level information and merging all
# processed samples into one single data object for the downstream of the analysis.
