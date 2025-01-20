##########################
## Loading packages and configs ##
##########################
suppressPackageStartupMessages({
  library(tidyverse)
  library(conflicted)
  library(qs)
  library(Seurat)
})
################################################################################
## Starting with the loading of the end product of "03main_nuc_frac_calcul.R" ##
################################################################################
big_sc <-
  qread(file = "data/non_qc-trimmed_nf-based_cell_status_added_big_sc.qs", nthreads = 10)

# We are not going to utilize "spliced" and "unspliced" transcript assays
# Reducing the file size of the object, with chopping the redundant bits

big_sc <-
  DietSeurat(
    object = big_sc,
    counts = T,
    scale.data = F,
    dimreducs = NULL,
    graphs = NULL,
    assays = c("RNA", "ADT")
  )
# During biological sample preparation, 4 emulsions were unintentionally handled to contain
# only T cells and only NK cells respectively. They were sequenced in separate emulsions.

# =====================================================================================================================================================
# table(big_sc$orig.ident)
#
# H21b       H21m       H22b       H22m       H23b       H23m       H24b       H24m       H25b       H25m       H26b       H26m       H27b       H27m
# 1578       7624        330      11467        397        554       2353       6081       2133       4683       4554       6922       8579       4431
# H28b       H28m       H29b       H29m       H30b       H30m    T22467b    T22467m T23491b_NK  T23491b_T T23491m_NK  T23491m_T    T23819b    T23819m
# 5679       9408      18914      11934       6510       4673       9532      12328        739      13931       1843       6960       1276       8218
# T23833b    T23833m    T24832b    T24832m    T25370b    T25370m    T25446b    T25446m    T26029b    T26029m    T26305b    T26305m    T26625b    T26625m
# 11008       5705       6397       3918       2982       9277       3336       3713      10940       9865      10350       3276       7987       5135
# =====================================================================================================================================================

# Using donor name only and embedding it into id of those samples
big_sc$id <-
  ifelse(
    big_sc$orig.ident == "T23491m_NK" |
      big_sc$orig.ident == "T23491m_T" |
      big_sc$orig.ident == "T23491b_NK" |
      big_sc$orig.ident == "T23491b_T",
    "T23491",
    str_sub(big_sc$orig.ident, 1, nchar(big_sc$orig.ident) - 1)
  )

# =====================================================================================================================================================
# > table(big_sc$id)
#
# H21    H22    H23    H24    H25    H26    H27    H28    H29    H30 T22467 T23491 T23819 T23833 T24832 T25370 T25446 T26029 T26305 T26625
# 9202  11797    951   8434   6816  11476  13010  15087  30848  11183  21860  23473   9494  16713  10315  12259   7049  20805  13626  13122
# =====================================================================================================================================================

##### Before QC-based trim plots ####
## #create_directory
before_trim_path <- "../pub_ready/docs/before_trim_QC_plots/"
dir.create(path = before_trim_path)
## nCountRNA faceted for 3 types of determined cell status per donor ##
p <-
  ggplot(big_sc@meta.data,
         aes(cell_status, nCount_RNA, fill = scDblFinder.class)) +
  geom_boxplot(outlier.size = 0) +
  facet_grid(cell_status ~ id, scales = "free", shrink = FALSE) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
# Saving this plot
ggsave(
  filename = "nCountRNA_faceted_cell_status.png",
  p,
  device = "png",
  dpi = 150,
  path = before_trim_path,
  height = 12,
  width = 16
)

## nFeatureRNA faceted for 3 types of determined cell status per donor ##
p <-
  ggplot(big_sc@meta.data,
         aes(cell_status, nFeature_RNA, fill = scDblFinder.class)) +
  geom_boxplot(outlier.size = 0) +
  facet_grid(cell_status ~ id, scales = "free") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
# Saving this plot
ggsave(
  filename = "nFeatureRNA_faceted_cell_status.png",
  p,
  device = "png",
  dpi = 150,
  path = before_trim_path,
  height = 12,
  width = 16
)

## percent.mt faceted for 3 types of determined cell status per donor ##
p <-
  ggplot(big_sc@meta.data,
         aes(cell_status, percent.mt, fill = scDblFinder.class)) +
  geom_boxplot(outlier.size = 0) +
  facet_grid(cell_status ~ id, scales = "free") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())#
# Saving this plot
ggsave(
  filename = "mito_percent_faceted_cell_status.png",
  p,
  device = "png",
  dpi = 150,
  path = before_trim_path,
  height = 12,
  width = 16
)


##### This is the decided QC trim criteria ####
meta.data <- big_sc@meta.data %>%
  dplyr::filter(cell_status == "cell") %>% #remove damaged cells and empty cells
  dplyr::filter(scDblFinder.class == "singlet") %>% #remove doublets
  dplyr::filter(nCount_RNA >= 500) %>% # remove bad quality cells
  dplyr::filter(nFeature_RNA >= 500) %>% # detected transcript types
  dplyr::filter(nCount_RNA <= 12000) %>% # detected transcript counts
  dplyr::filter(percent.mt <= 10) %>% #percentage of mitochondrial gene transcripts
  dplyr::filter(predicted.celltype.l1.score >= 0.90) # Prediction confidence of main immune cell type

## #create_directory
after_trim_path <- "../pub_ready/docs/after_trim_QC_plots/"
dir.create(path = after_trim_path)
## nCountRNA faceted for 3 types of determined cell status per donor ##
p <-
  ggplot(meta.data,
         aes(cell_status, nCount_RNA, fill = scDblFinder.class)) +
  geom_boxplot(outlier.size = 0) +
  facet_grid(cell_status ~ id, scales = "free", shrink = FALSE) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
# Saving this plot
ggsave(
  filename = "nCountRNA_faceted_cell_status.png",
  p,
  device = "png",
  dpi = 150,
  path = after_trim_path,
  height = 12,
  width = 16
)

## nFeatureRNA faceted for 3 types of determined cell status per donor ##
p <-
  ggplot(meta.data,
         aes(cell_status, nFeature_RNA, fill = scDblFinder.class)) +
  geom_boxplot(outlier.size = 0) +
  facet_grid(cell_status ~ id, scales = "free") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
# Saving this plot
ggsave(
  filename = "nFeatureRNA_faceted_cell_status.png",
  p,
  device = "png",
  dpi = 150,
  path = after_trim_path,
  height = 12,
  width = 16
)

## percent.mt faceted for 3 types of determined cell status per donor ##
p <-
  ggplot(meta.data,
         aes(cell_status, percent.mt, fill = scDblFinder.class)) +
  geom_boxplot(outlier.size = 0) +
  facet_grid(cell_status ~ id, scales = "free") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())#
# Saving this plot
ggsave(
  filename = "mito_percent_faceted_cell_status.png",
  p,
  device = "png",
  dpi = 150,
  path = after_trim_path,
  height = 12,
  width = 16
)


trim_sc <-
  subset(
    big_sc,
    cell_status == "cell" & #remove damaged cells and empty cells
      scDblFinder.class == "singlet" &
      #remove doublets
      nCount_RNA >= 500 &
      #remove bad quality cells
      nFeature_RNA >= 500 &
      #detected transcript types
      nCount_RNA <= 12000 &
      #detected total transcript counts
      percent.mt <= 10 &
      #percentage of mitochondrial gene transcripts
      predicted.celltype.l1.score >= 0.90
  ) #Prediction confidence of main immune cell type

trim_sc
# An object of class Seurat
# 36615 features across 109690 samples within 2 assays
# Active assay: RNA (36601 features, 0 variable features)
# 1 other assay present: ADT
table(trim_sc$orig.ident)
# Good-quality of the cells per donor emulsion

# > table(trim_sc$orig.ident)
#
# H21b       H21m       H22b       H22m       H23b       H23m       H24b       H24m       H25b       H25m       H26b       H26m       H27b       H27m
# 918       4546        178       6612         67        293        529       2306        985       1621       2470       4233       5551       2770
# H28b       H28m       H29b       H29m       H30b       H30m    T22467b    T22467m T23491b_NK  T23491b_T T23491m_NK  T23491m_T    T23819b    T23819m
# 2692       4031       2950       3916       2177       1811       5133       4549        350       6012        896       3071        482       2782
# T23833b    T23833m    T24832b    T24832m    T25370b    T25370m    T25446b    T25446m    T26029b    T26029m    T26305b    T26305m    T26625b    T26625m
# 4179       1871       2052       1873       1607       3766       1720       1840       5302       4220       2938        522       2114       1755


qsave(trim_sc,file = "data/qc_trimmed_t_plus_nk.qs")
