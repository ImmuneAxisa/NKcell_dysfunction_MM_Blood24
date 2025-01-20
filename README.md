## Summary

This repository contains the scripts yielding to the processed scRNAseq data of NK cells in Multiple Myeloma patients.
The study utilizing this dataset is published in   the results in our [article](https://doi.org/10.1182/blood.2023023529).

The processed data of 8 NK cell clusters determined in the MM patients and healthy donors were stored in [Zenodo](https://doi.org/10.5281/zenodo.13359147)
![Screenshot from 2024-08-22 13-56-34](https://github.com/user-attachments/assets/e01e16c4-1030-4a1e-858a-edab776095a1)

Below structure of scripts are the step by step data curation and then investigation of NK cells in the downstream analysis.

## Analysis diagram

![AnalysisDiagram](https://github.com/user-attachments/assets/28c07e18-66f0-408e-9213-465b946abfba)


## Structural layout

only code is stored on github but full structure with intermediate files can be retrieved on the [Zenodo repository](https://doi.org/10.5281/zenodo.13359147)

```{r eval=FALSE, include=TRUE}
NKcell_dysfunction_MM_Blood24
├── LICENSE
├── NKcell_dysfunction_MM_Blood24.Rproj
├── README.md #this file
├── index.html # rendered github pages website index
├── index.md #github pages website index
├── data
│   └── List Pathways_clusters.xlsx #curated genesets from publications
├── doc
│   ├── NKanalysis.Rmd # initial exploratory analysis
│   ├── NKanalysis.html # rendered initial exploratory analysis
│   ├── RNAseq_figs.Rmd # advanced analysis and publication figures construction
│   └── RNAseq_figs.html # rendered advanced analysis and publication figures construction
├── figs
│   ├── manuscript # pdf panels for manuscript figure construction
│   │   ├── MELD_UMAP.pdf
│   │   ├── MELD_violins.pdf
│   │   ├── custom_genesets_hm.pdf
│   │   ├── disease_hm_kegg.pdf
│   │   ├── disease_leadingEdge_hm.pdf
│   │   ├── frequency.pdf
│   │   ├── gsea_pair_results.pdf
│   │   ├── gsea_pair_results_reactome.pdf
│   │   ├── gsea_results_2.xlsx
│   │   ├── leadingEdge_hm.pdf
│   │   ├── main_umaps.pdf
│   │   ├── patient_hm_kegg.pdf
│   │   ├── pseudotime_results.pdf
│   │   ├── pt_genes.pdf
│   │   ├── slingshot_C8_start.pdf
│   │   ├── slingshot_bmmc.pdf
│   │   └── slingshot_pbmc.pdf
│   └── markers
│       └── NK # folder containing feature plots of relevant genes in png format
├── file_tree.txt
├── results
│   └── NK_analysis
│       ├── clustering_typing # tabular data of analysis outputs
│       │   ├── clusters_markers_table_harmony_NK.qs
│       │   ├── clusters_markers_table_harmony_NK_filtered.qs
│       │   ├── clusters_markers_table_harmony_NK_paired.qs
│       │   ├── clusters_markers_table_harmony_NK_paired_filtered.csv
│       │   ├── clusters_markers_table_harmony_NK_paired_filtered.qs
│       │   ├── clusters_markers_table_harmony_NK_paired_filtered_summarized.csv
│       │   ├── clusters_markers_table_paired_from_ruchan_function.qs
│       │   ├── custom_NK_genesets.qs
│       │   ├── custom_NK_genesets_cleaned.qs
│       │   └── gsea_results.qs
│       ├── phate_meld_output
│       │   ├── HDvsMM_BMMC_NK_MELD_LLH.csv
│       │   ├── HDvsMM_PBMC_NK_MELD_LLH.csv
│       │   ├── HDvsMM_PBMC_NK_MELD_densities.csv
│       │   ├── PBMCvsBMMC_HD_NK_MELD_LLH.csv
│       │   └── PBMCvsBMMC_MM_NK_MELD_LLH.csv
│       ├── seurat_objects
│       │   ├── analysisNK_Ruchan.qs
│       │   └── seurat_diet_analysisNK_Ruchan.rds
│       └── trajectory_pseudotime
│           ├── slingshot_pto.qs # main analysis
│           ├── slingshot_pto_c5_bmmc.qs # analysis split by tissue
│           ├── slingshot_pto_c5_pbmc.qs # analysis split by tissue
│           └── slingshot_pto_c8.qs # cluster 8 start
└── src
    ├── MELD # see README for details
    │   ├── MELD_v2.1.py
    │   ├── PHATE_v2.py
    │   ├── README.html
    │   ├── README.md
    │   ├── phate_meld_sbatcher_w_args_v2.sh
    │   └── phate_meld_wrapper_v2_NK_analysis.R
    ├── clustering_typing 
    │   └── C01_cluster_markers_NK.R # get most relevant markers for clusters
    ├── data_preparation_and_QC
    │   ├── 01_individually_processing_per_emulsion.R # Reading Cellranger / Velocyto outputs, then annotating cell types per emulsion to ref PBMC in a loop
    │   ├── 02_merging_processed_emulsions.R #Merging each emulsions' data into one data object
    │   ├── 03a_nuc_frac_calcul.R #Calculating nuclear RNA fraction of cells
    │   ├── 03b_nuclear_RNA_fraction_metric.py # Helper python script of nuclear RNA fraction
    │   ├── 04_QC_check_and_trim.R # QC plots and trimming based-on QC
    │   ├── 05a_preprocessing.R # Normalization,ScaleData,Dimensional Reduction,PC-correction etc.
    │   ├── 05b_plots_leading_to_nk_subsetting.R # Intermediate step of assessing inferred cell types
    │   ├── 06a_nk_subsetting.R # Separating NK cells from the rest of the cell types 
    │   ├── 06b_removal_of_prolif_nk_portion.R # Keeping non-proliferating, more stable NK cells for downstream analysis
    │   ├── combin_gsea_function.R # Multiple parameter assessment oriented GSEA script
    │   ├── combin_presto_function.R # Helper function of combinatorial Presto running to rank all genes for comparison groups (ClusterXvsY)
    │   └── multi_criteria_filt_GSEA.R # Helper function to run fGSEA for comparison pairs and keep them in list of lists format
    ├── data_retrieval # get objects formatted to use within project
    │   ├── diet_seurat_for_MELD.R
    │   └── get_robjects_from_server.R
    ├── helpers # plotting helpers
    │   └── slingshot_layout_helpers.R
    └── trajectory # slingshot scripts
        ├── slingshot_c8_start.R
        └── slingshot_sep_organs.R
```

## Analysis medium

All of the analyses of this study except *MELD* were performed with the local machine (workstation) containing Intel® Core™ i9-9900K CPU and 128 GB RAM, running on Linux (Ubuntu 22.04.2 LTS) operation system.
