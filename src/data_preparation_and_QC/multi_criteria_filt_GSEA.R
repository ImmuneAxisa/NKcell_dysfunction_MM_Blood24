######## suppressPackageStartupMessages() is a function that suppresses the messages that are printed when a package is loaded.
# library() is a function that loads a package.
suppressPackageStartupMessages({
        library(qs) # loads the qs package
        library(Seurat) # loads the Seurat package
        library(presto) # loads the presto package
        library(msigdbr) # loads the msigdbr package
        library(fgsea) # loads the fgsea package
        library(dplyr) # loads the dplyr package
        library(ggplot2) # loads the ggplot2 package
        library(tidyverse) # loads the tidyverse package
        library(scran) # loads the scran package
        library(here) # loads the here package
})

print(here()) # prints the current working directory

#### Reading the sc_obj ####
sc_obj <-qread(here("data", "nonprolif_NKs_multiple_res_clusters_20dec22.qs"))

#### Running Pairwise Presto by all combination pairs ##
source(here("code/", "combin_presto_function.R"))
# Calculating the all pair-wise presto results of the clusters
paired_results <- combin_presto(sc_obj = sc_obj) # 62.065 sec elapsed
gc()
# **************************************************************************

##### Calling  FGSEA Custom Function here ##
source(here("code/", "combin_gsea_function.R"))

gene_counts <- rowSums(sc_obj@assays$RNA@counts)
length(gene_counts[gene_counts < 1]) # 6528 genes has less than 1 detection
gene_detection <- rowSums(sc_obj@assays$RNA@counts > 0) / ncol(sc_obj) *100 #calculating the per gene detection rate
qplot(x = gene_detection) +stat_ecdf(aes(after_stat(y) * 100))+ scale_x_continuous(trans = "log10") + annotation_logticks()+ geom_vline(aes(xintercept = 2,colour ="red"))+NoLegend()

length(gene_detection[gene_detection > 1])

subset_sc_det1 <- sc_obj[gene_detection > 1]

subset_det1_paired_res <- combin_presto(sc_obj = subset_sc_det1)

########### Running for all combins

det_thr_1_auc_all <- run_combin_GSEA(paired_results = subset_det1_paired_res,rank_how = "auc",genesets_db = genesets_db)
qsave(det_thr_1_auc_all,file = "detection_1_thresholded_auc_all_combin_gsea_param1_res.qs")

no_thr_auc_all <- run_combin_GSEA(paired_results = paired_results,rank_how = "auc",genesets_db = genesets_db)
qsave(no_thr_auc_all,file = "non_thresholded_auc_all_combin_gsea_param1_res.qs")


det_thr_1_logfc_all <- run_combin_GSEA(paired_results = subset_det1_paired_res,rank_how = "logfc",genesets_db = genesets_db)
qsave(det_thr_1_logfc_all,file = "detection_1_thresholded_logfc_all_combin_gsea_param1_res.qs")

no_thr_logfc_all <- run_combin_GSEA(paired_results = paired_results,rank_how = "logfc",genesets_db = genesets_db)
qsave(no_thr_logfc_all,file = "non_thresholded_logfc_all_combin_gsea_param1_res.qs")

####################
####################

no_thr_auc_all <- clean_up_gsea(no_thr_auc_all,ranking_metric = "auc",thresholding = "no_thr")
no_thr_logfc_all <- clean_up_gsea(no_thr_logfc_all,ranking_metric = "logfc",thresholding = "no_thr")

det_thr_1_auc_all <- clean_up_gsea(det_thr_1_auc_all,ranking_metric = "auc",thresholding = "det_1") # detection rate of expressed cells 1 in each cell
det_thr_1_logfc_all <- clean_up_gsea(det_thr_1_logfc_all,ranking_metric = "logfc",thresholding = "det_1") # detection rate of expressed cells 1 in each cell



## do the  rbind of all iteration


res <- rbind(no_thr_auc_all,det_thr_1_auc_all,no_thr_logfc_all,det_thr_1_logfc_all)


## calculate correlation across iterations

cor_combin <- res %>%
        # mutate(ranking = "auc",
        #        threshold = "no_thr") %>%
        unite("grouping", ranking, threshold) %>%
        pivot_wider(names_from = grouping, values_from = NES) %>%
        group_by(database, pair) %>%
        nest() %>%
        mutate(correlation = purrr::map(data, ~ .x %>%
                                                dplyr::select(-pathway) %>%
                                                correlate() %>%
                                                stretch(na.rm = TRUE, remove.dups =FALSE))) %>%
        dplyr::select(-data) %>%
        unnest(cols = c(correlation)) %>%
        dplyr::filter(x != y)


ggplot(data = cor_combin) +
        geom_point(mapping = aes(x = r, y = r, color = database,
                                 shape = x))

ggplot(data = cor_combin) + geom_point(aes(x = pair, y=x ,color = r)) + scale_colour_viridis_b() + facet_grid( database ~ y) + theme_dark(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot(data = cor_combin, mapping = aes(x = x, y = r)) + geom_boxplot(show.legend = T) +
        scale_colour_viridis_b() + facet_grid(database ~ y) + #theme_minimal(base_size = 14) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
