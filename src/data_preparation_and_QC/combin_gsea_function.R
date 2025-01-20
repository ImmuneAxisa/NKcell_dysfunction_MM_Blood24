suppressPackageStartupMessages({
 library(qs)
 library(Seurat)
 library(presto)
 library(msigdbr)
 library(fgsea)
 library(dplyr)
 library(ggplot2)
 library(tidyverse)
 library(scran)
 library(here)
 library(tictoc)
 library(purrr)
})



#### check/make the frozen GSEA database ####
# check if the file exists
ifelse (test = file.exists(here("data/GSEA_DB","static_snapshot.qs")),
        yes = {
                # load the file
                genesets_db <- qread(file = here("data/GSEA_DB", "static_snapshot.qs"))
                # print the date of the file
                print("Loading the selected GSEA sets from January 2023")

        },
        no = {
                # print that the file is empty
                print("static_snapshot file is empty")
                tictoc::tic()
                # if the file does not exist, create the file
                # GO Biological process
                bio_pro <-
                        msigdbr(
                                species = "Homo sapiens",
                                category = "C5",
                                subcategory = "GO:BP"
                        )
                bio_pro <- bio_pro %>% split(., x = .$gene_symbol, f = .$gs_name)
                #reactome
                react <-
                        msigdbr(
                                species = "Homo sapiens",
                                category = "C2",
                                subcategory = "REACTOME"
                        )
                react <- react %>% split(., x = .$gene_symbol, f = .$gs_name)
                #Transcription factor targets
                tft <-
                        msigdbr(
                                species = "Homo sapiens",
                                category = "C3",
                                subcategory = "TFT:GTRD"
                        )
                tft <- tft %>% split(., x = .$gene_symbol, f = .$gs_name)
                #kegg
                kegg <-
                        msigdbr(
                                species = "Homo sapiens",
                                category = "C2",
                                subcategory = "KEGG"
                        )
                kegg <- kegg %>% split(., x = .$gene_symbol, f = .$gs_name)
                # create a list of the genesets
                genesets_db <-
                        list(kegg, bio_pro, react, tft) %>% set_names("kegg", "bio_pro", "react", "tft")
                # genesets_db <- genesets_db %>% map( ~ split(., .x$gene_symbol, f = .x$gs_name))
                # remove the genesets from the environment
                rm(kegg); rm(bio_pro); rm(react); rm(tft)
                # save the genesets to the file
                qsave(x = genesets_db,
                      file = here("data/GSEA_DB", "static_snapshot.qs"))
                # print that the file is saved
                print(paste0(
                        "static_snapshot of GSEA DB is saved into",
                        here("data/GSEA_DB", "static_snapshot.qs")
                ))
                # load the file
                genesets_db <- qread(file = here("data/GSEA_DB", "static_snapshot.qs"))
                # print the date of the file
                print("and then loaded the selected GSEA sets from January 2023")
                tictoc::toc()
        })
#####
run_combin_GSEA <- function(paired_results, rank_how, genesets_db) {
        set.seed(2206)
        results <- list()
        for (i in 1:length(paired_results)) {
                if (rank_how == "log_fc" || rank_how == "logfc" || rank_how == "logFC" || rank_how == "log_FC" || rank_how == "log_Fc") {
                        ranked_genes <- paired_results[[i]] %>% dplyr::select(feature, logFC)
                } else if (rank_how == "auc" || rank_how == "AuC" || rank_how == "AUC") {
                        ranked_genes <- paired_results[[i]] %>% dplyr::select(feature, auc)
                        ranked_genes$auc <- ranked_genes$auc-0.5 #centering AuC to Zero
                } else {
                        stop("Non-defined ranking, sorry !")
                }
                ranked_genes <- deframe(ranked_genes)
                print(ranked_genes)
                fgsea_results <- lapply(genesets_db, function(y) {
                        fgseaRes <- fgsea(pathways = y,
                                stats = ranked_genes,
                                minSize = 15,
                                maxSize = 500,
                                gseaParam = 1
                        )

                        return(fgseaRes)

                })


                #
                results[[i]] <- fgsea_results


                # attaching ranked_genes for plotting purposes of fgsea
                results[[i]]$ranks <- ranked_genes
                results[[i]]$pair <- unique(paired_results[[i]]$pairs)
        }
        return(results)
        }
#####
#####

#### Tidy clean up function  for nested
#### list structure of fgsea result embedddings

library(corrr)

# res <- qs::qread("/Volumes/CRCT13/Ruchan/01_NK_project/01_human/8_clusters_nonproliferating_nk/Enrichment_results/detection_1_thresholded_all_combin_gsea_res.qs")


clean_up_gsea <- function(res, ranking_metric, thresholding) {
        pair_names <- lapply(res, function(it)
                it[["pair"]])

        pair_names <- do.call(c, pair_names)


        names(res) <- pair_names


        res <- lapply(pair_names, function(pair_it) {
                tmp_res <- res[[pair_it]]

                tmp_res[["pair"]] <- NULL

                new_names <- paste0(names(tmp_res), "-", pair_it)

                names(tmp_res) <- new_names

                return(tmp_res)
        })


        res <- do.call(c, res)

        res <- enframe(res, "iteration", "results")

# Re-organizing an grouping per pairs
# in this task, *NES*: "Normalized Enrichment Score"
# is the interesting bit of comparison for us.
        res %>%
                separate(iteration, c("database", "pair"), sep = "-") %>%
                dplyr::filter(database != "ranks") %>%
                mutate(results = purrr::map(results, as.data.frame)) %>%
                unnest(cols = c(results)) %>%
                dplyr::select(database, pair, pathway, NES) %>%
                mutate(ranking = ranking_metric,
                       threshold = thresholding)


}
