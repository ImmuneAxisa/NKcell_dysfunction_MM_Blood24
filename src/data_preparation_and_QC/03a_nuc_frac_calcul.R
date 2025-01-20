if(!require(DropletQC)){
 install.packages("DropletQC")
 library(DropletQC)
}
if(!require(anndata)){
 install.packages("anndata")
 library(anndata)
}
library(reticulate)
library(qs)
library(dplyr)
setwd(dir = "/mnt/depot/pub_ready/")
py_run_file(file = "03b_nuclear_RNA_fraction_metric.py",convert = TRUE,local = FALSE)
tmp <- py$adata
cell_barcodes <- tidyr::as_tibble(rownames(tmp))
cell_barcodes <- tidyr::as_tibble(substr(cell_barcodes$value,1, nchar(cell_barcodes$value)-1))
nf <- py$nuclear_fraction
nf_ann <- cbind(cell_barcodes,nf)
names(nf_ann) <- c("cell_barcode","nuc_frac")
qs::qsave(nf_ann,"data/nf_ann.qs")

big_sc <- qs::qread(file = "data/velo_embedded_big_sc_whole_pbmc_annotated.qs")
metdat <- big_sc@meta.data
metdat <- rownames_to_column(.data = metdat,var = "cell_barcode")

qs::qsave(metdat,file = "data/non_qc_trimmed_nf_added_big_sc_metadata.qs")


#Read the nf_added metadata here
metdat <- qread(file = "data/non_qc_trimmed_nf_added_big_sc_metadata.qs")
# Calculate the nuclear fraction using the spliced and unspliced matrices
# tmp_ed <- identify_empty_drops(nf_umi =as.data.frame(cbind(metdat$nuc_frac,metdat$nCount_RNA)), include_plot = T,
#                                plot_path = "docs/nuc_frac_plots/",pdf_png = "png",plot_name = "identified_empty_droplets.png",
#                                plot_width = 16,plot_height = 12)
#
#
# Customization of the internal function of the identify_empty_drops from DropletQC package
# https://github.com/powellgenomicslab/DropletQC/blob/main/R/identify_empty_drops.R
custom_identify_empty_drops <-
  function(nf_umi,
           nf_rescue = 0.05,
           umi_rescue = 1000,
           include_plot = FALSE,
           plot_name = NULL,
           plot_path = NULL,
           plot_width=18,
           plot_height=13,
           pdf_png = "png"
  ) {
    library(patchwork)
    ## Check and parse arguments
    if (any(class(nf_umi) == "data.frame")) {

      # Assume nuclear fraction is in the first column
      nf <- unlist(nf_umi[, 1], use.names = FALSE)
      # Assume UMI counts are in the second column
      umi <- unlist(nf_umi[, 2], use.names = FALSE)

      # Check values are reasonable
      if(any(c(max(nf)>1, min(nf)<0))){
        warning(paste0("The nuclear fraction values provided in the first column of 'nf_umi' should be between 0 and 1, but values outside this range were identified : minimum = ",min(nf),", maximum = ",max(nf)), call.=FALSE)
      }

      if(!all(umi == floor(umi))){
        non_integer_examples <- which(umi != floor(umi))
        if(length(non_integer_examples)>5){
          non_integer_examples <- non_integer_examples[1:5]
        }
        non_integer_examples <- paste(umi[non_integer_examples], collapse = ",")
        warning(paste0("Non-integer values detected in the second column of 'nf_umi' (e.g. ",non_integer_examples,") where umi counts were expected"), call.=FALSE)
      }

      if(max(umi)<100){
        umi
        warning(paste0("The total umi counts provided in the second column of 'nf_umi' appear to be quite low (max = ",max(umi),"), are these the total UMI counts per cell?"), call.=FALSE)
      }

    } else {
      stop(paste0("A data frame should be supplied to the nf_umi argument, but an object of class ",paste(class(nf_umi), collapse = "/")," was provided"), call.=FALSE)
    }

    # Density estimation (automatically chosen bandwidth)
    kdde_0 <- ks::kdde(x = nf, deriv.order = 0)
    kdde_0 <- data.frame(estimate = kdde_0[["estimate"]],
                         eval.points = kdde_0[["eval.points"]])

    # Density derivative estimation (automatically chosen bandwidth, but different
    # from kdde_0!)
    kdde_1 <- ks::kdde(x = nf, deriv.order = 1)
    kdde_1 <- data.frame(estimate = kdde_1[["estimate"]],
                         eval.points = kdde_1[["eval.points"]])

    # Find point to place cut-off between empty droplets and cells
    gradient_sign <- rle(kdde_1[["estimate"]]>0)
    nf_cutoff <- kdde_1[["eval.points"]][sum(gradient_sign[["lengths"]][1:2])]

    ## Check if there is more than one peak
    if(length(gradient_sign$values)<4){
      warning(paste0("Could not detect more than one peak in the nuclear fraction distribution. There may not be any empty droplets present. We suggest visualising the density estimation (include_plot=TRUE)."), call.=FALSE)
    }
    library(ggplot2)
    # Label cells
    nf_umi$cell_status <- "cell"
    nf_umi$cell_status[nf < nf_cutoff] <- "empty_droplet"

    # Rescue miscalled cells
    nf_umi$cell_status[nf > nf_rescue &  umi > umi_rescue] <- "cell"

    # Plots
    if (include_plot) {

      p1 <- ggplot2::ggplot(kdde_0, ggplot2::aes(x = eval.points, y = estimate)) +
        ggplot2::geom_line() +
        ggplot2::ggtitle("Density estimation") +
        ggplot2::xlab("Nuclear fraction") +
        ggplot2::ylab("Density function") +
        ggplot2::geom_vline(xintercept = nf_cutoff,
                            col = "dodgerblue", linetype = "dashed")


      p2 <- ggplot2::ggplot(kdde_1, ggplot2::aes(x = eval.points, y = estimate)) +
        ggplot2::geom_line() +
        ggplot2::ggtitle("Density derivative estimation") +
        ggplot2::xlab("Nuclear fraction") +
        ggplot2::ylab("Density derivative function") +
        ggplot2::geom_vline(xintercept = nf_cutoff,
                            col = "dodgerblue", linetype = "dashed") +
        ggplot2::geom_hline(yintercept = 0, col = "grey")

      p3_data <- data.frame(nf = nf,
                            umi = umi,
                            cell_status = nf_umi$cell_status)
      p3 <- ggplot2::ggplot(p3_data,
                            ggplot2::aes(
                              x = nf,
                              y = log10(umi),
                              colour = cell_status
                            )) +
        ggplot2::geom_point(size = 0.3,alpha = 0.3) +  # Adjust point size as needed
        ggplot2::xlim(0, 1) +
        ggplot2::xlab("Nuclear fraction") +
        ggplot2::ylab("log10(UMI count)") +
        ggplot2::theme(legend.title = ggplot2::element_blank()) +
        ggplot2::guides(colour = ggplot2::guide_legend(title = "Cell Status")) +
        ggplot2::geom_density_2d(aes(color = cell_status), alpha = 0.1) + # Add density embeddings
        ggplot2::scale_color_manual(values = c("cell" = "#4575b4", "empty_droplet" = "#565656"))


      p3
      # Add three dashed lines to indicate threshold + rescue area
      l1 <- data.frame(x_start = ifelse(umi_rescue<min(umi) & nf_rescue<nf_cutoff, max(min(nf),nf_rescue), nf_cutoff),
                       x_end = ifelse(umi_rescue<min(umi) & nf_rescue<nf_cutoff, max(min(nf),nf_rescue), nf_cutoff),
                       y_start = log10(min(umi)),
                       y_end = max(log10(umi_rescue), log10(min(umi))))

      l2 <- data.frame(x_start = ifelse(umi_rescue<min(umi) & nf_rescue<nf_cutoff, max(min(nf),nf_rescue), nf_cutoff),
                       x_end = min(max(min(nf),nf_rescue), nf_cutoff),
                       y_start = max(log10(umi_rescue), log10(min(umi))),
                       y_end = max(log10(umi_rescue), log10(min(umi))))

      l3 <- data.frame(x_start = min(max(min(nf),nf_rescue), nf_cutoff),
                       x_end = min(max(min(nf),nf_rescue), nf_cutoff),
                       y_start = max(log10(umi_rescue), log10(min(umi))),
                       y_end = max(log10(umi)))

      p3 <- p3 + ggplot2::geom_segment(ggplot2::aes(
        x = l1$x_start,
        y = l1$y_start,
        xend = l1$x_end,
        yend = l1$y_end),
        colour = "dodgerblue",
        linetype = "dashed",
        data = p3_data)

      p3 <- p3 + ggplot2::geom_segment(ggplot2::aes(
        x = l2$x_start,
        y = l2$y_start,
        xend = l2$x_end,
        yend = l2$y_end),
        colour = "dodgerblue",
        linetype = "dashed",
        data = p3_data)

      p3 <- p3 + ggplot2::geom_segment(ggplot2::aes(
        x = l3$x_start,
        y = l3$y_start,
        xend = l3$x_end,
        yend = l3$y_end),
        colour = "dodgerblue",
        linetype = "dashed",
        data = p3_data)


      # Combine plots using patchwork
      combined_plots <- (p1 | p2 | p3)

      # Save or print the combined plot
      if (is.null(plot_path)) {
        print(combined_plots)
      } else {
        ggsave(filename = plot_name,
               plot = combined_plots,
               device = pdf_png,
               path = plot_path,
               width = plot_width,
               height = plot_height,
               units = "cm")
}
    }

    return(nf_umi)
  }


tmp_ed <- custom_identify_empty_drops(nf_umi =as.data.frame(cbind(metdat$nuc_frac,metdat$nCount_RNA)), include_plot = T,
                                      plot_path = "docs/nuc_frac_plots/",pdf_png = "png",plot_name = "identified_empty_droplets.png",
                                      plot_width = 30,plot_height = 6)



## Nuclear RNA fraction differs per cell type
##
####### Customization injection to "identify_damaged_cells" function of DropletQC package
# https://github.com/powellgenomicslab/DropletQC/blob/main/R/identify_damaged_cells.R

# Generating Gaussian mixture model with a maximum of two components and fitting to the umi counts and nuclear fraction scores per cell type.
# The parameters of the model are estimated using expectation maximisation (EM) with the mclust package.
# The best model is selected using the Bayesian Information Criterion.
# The two populations (cells and damaged cells) are assumed to have equal variance.
# Finally, labeling determined damaged cells and "whole" cells in to "cell_status".
# ***************
# tmp_dc <-
#  identify_damaged_cells(
#   nf_umi_ed_ct = as_data_frame(cbind(tmp_ed, (
#    big_sc$predicted.celltype.l2
#   ))),
#   nf_sep = 0.05,
#   umi_sep_perc = 15,
#   output_plots = TRUE
#  )
# **********************************
assess_EM <- function(em, umi_thresh, nf_thresh){

  # Separately for each cell type, assign a barcode as "cell" or "damaged_cell"
  # based on the EM results using the following sequential procedure:

  # 1. If the EM model selected contained only one distribution, score all
  # barcodes as "cell".
  check_1 <- em$G==2

  # 2. If two distributions were fit, check the distribution with the higher
  # nuclear_fraction mean also has a lower umi mean. If this criteria is
  # satisfied, we consider the population with the lower umi count and higher
  # nuclear_fraction scores the damaged_cell population and move on to step 3.
  # If this is not the case we score all barcodes as "cell".
  if(check_1){
    nf_means <- em$parameters$mean["nf",]
    umi_means <- em$parameters$mean["umi",]
    check_2 <- umi_means[which.max(nf_means)] < umi_means[which.min(nf_means)]
  } else {
    check_2 <- FALSE
  }

  # 3. If 2. was satisfied, we check the damaged_cell nuclear_fraction
  # distribution mean is at least nf_thresh greater than the cell mean. Also
  # check the damaged_cell umi distribution mean is at least umi_thresh
  # percent lower than the cell umi distribution mean. If these two criteria
  # are satisfied we assign the damaged cells, otherwise all barcodes are
  # labelled "cell".
  if(check_2){
    # Check nuclear fraction threshold is satisfied
    nf_check <- nf_means[which.max(nf_means)] - nf_means[which.min(nf_means)] > nf_thresh
    # Check umi threshold is satisfied
    damaged_cell_umi <- 10^umi_means[which.max(nf_means)]
    cell_umi <- 10^umi_means[which.min(nf_means)]
    umi_check <- damaged_cell_umi < (cell_umi - cell_umi*(umi_thresh/100))

    if(all(nf_check, umi_check)){ check_3 <- TRUE } else { check_3 <- FALSE }
  } else {
    check_3 <- FALSE
  }

  # Assign barcodes
  em_classification <- data.frame(nf = em$data[,"nf"],
                                  umi = em$data[,"umi"],
                                  classification= em$classification)
  row.names(em_classification) <- row.names(em$data)

  if(all(check_1, check_2, check_3)){
    damaged_cells <- em_classification$classification==which.max(em$parameters$mean["nf",])
    em_classification$classification[damaged_cells] <- "damaged_cell"
    em_classification$classification[!damaged_cells] <- "cell"
  } else {
    em_classification$classification <- "cell"
  }

  return(em_classification)

}

#' Plot EM results
#'
#' @description This is a helper function called internally by the
#'   `identify_damaged_cells` function (if plots are requested). It's not
#'   intended for general use.
#'
#' @param emMclust Mclust, result of running EM on the log10(UMI counts) and the
#'   nuclear fraction
#' @param em_classified data frame, with three columns; log10(UMI counts),
#'   nuclear fraction and a column defining each cell as a "cell" or
#'   "damaged_cell" - this should be the output from running `assess_EM` on the
#'   model results
#' @param input_cell_type character, the name of the cell type
#'
#' @return ggarrange, three ggplots combined using `ggpubr::ggarrange`
#'
#' @keywords internal
plot_EM <- function(em, em_c, input_cell_type, umi_thresh, nf_thresh){

  # Get model parameters for plotting
  em_params <- list(nf_means = em$parameters$mean["nf",],
                    umi_means = em$parameters$mean["umi",],
                    nf_sd = rep(sqrt(em$parameters$variance$Sigma["nf",1]), em$G),
                    umi_sd = rep(sqrt(em$parameters$variance$Sigma["umi",2]), em$G),
                    amplitude = em$parameters$pro)
  # Add colours
  cell_colour <- "#4575b4"
  damaged_cell_colour <- "darkred"
  if(em$G==2){
    em_params$distribution_colour <- c(damaged_cell_colour, cell_colour)[order(em_params$nf_means, decreasing = TRUE)]
  } else {
    em_params$distribution_colour <- cell_colour
  }

  # Plot

  ### nf vs. umi ###
  p1 <- ggplot2::ggplot(em_c, ggplot2::aes(x=nf, y=umi, colour=classification)) +
    ggplot2::geom_point(size = 0.5,alpha = 0.3) +
    ggplot2::xlab("Nuclear fraction") +
    ggplot2::ylab("log10(UMI count)") +
    ggplot2::guides(colour = ggplot2::guide_legend(title = "Cell Status")) +
    ggplot2::scale_color_manual(values = c("cell" = "#4575b4", "damaged_cell" = "darkred"))

  ### nf ###
  p2 <- ggplot2::ggplot(em_c, ggplot2::aes(x = nf)) +
    ggplot2::geom_histogram(breaks = seq(0, 1, len = 51), colour = "black", fill = "white") +
    mapply(
      function(mean, sd, amplitude, n, binwidth, dist_colour) {
        ggplot2::stat_function(
          fun = function(x) {
            (stats::dnorm(x, mean = mean, sd = sd)) * n * binwidth * amplitude
          },
          size=1.5, geom="area", alpha=0.5, fill=dist_colour)
      },
      mean = em_params[["nf_means"]], #mean
      sd = em_params[["nf_sd"]], #standard deviation
      amplitude = em_params[["amplitude"]], #amplitude
      n = nrow(em_c), #sample size
      binwidth = 1/51, #binwidth used for histogram
      dist_colour = em_params[["distribution_colour"]]
    ) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab("Nuclear fraction")
  # add lines for damaged cell nuclear fraction mean (solid, red) and threshold (dashed, blue)
  if(em$G==2){
    p2 <- p2 + ggplot2::geom_vline(xintercept = em_params[["nf_means"]][which.max(em_params[["nf_means"]])], colour = damaged_cell_colour)
    p2 <- p2 + ggplot2::geom_vline(xintercept = em_params[["nf_means"]][which.min(em_params[["nf_means"]])] + nf_thresh, linetype="dashed", colour = cell_colour)
  }


  ### umi ###
  p3 <- ggplot2::ggplot(em_c, ggplot2::aes(x = umi)) +
    ggplot2::geom_histogram(breaks = seq(min(em_c$umi), max(em_c$umi), len = 51), colour = "black", fill = "white") +
    mapply(
      function(mean, sd, amplitude, n, binwidth, dist_colour) {
        ggplot2::stat_function(
          fun = function(x) {
            (stats::dnorm(x, mean = mean, sd = sd)) * n * binwidth * amplitude
          },
          size=1.5, geom="area", alpha=0.5, fill=dist_colour)
      },
      mean = em_params[["umi_means"]], #mean
      sd = em_params[["umi_sd"]], #standard deviation
      amplitude = em_params[["amplitude"]], #amplitude
      n = nrow(em_c), #sample size
      binwidth = c(max(em_c$umi) - min(em_c$umi))/51, #binwidth used for histogram
      dist_colour = em_params[["distribution_colour"]]
    ) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab("log10(UMI count)")
  # add lines for damaged cell umi mean (solid, red) and threshold (dashed, blue)
  if(em$G==2){
    # add line for damaged cell umi mean (solid, red)
    p3 <- p3 + ggplot2::geom_vline(xintercept = em_params[["umi_means"]][which.max(em_params[["nf_means"]])], colour = damaged_cell_colour)
    # add lines for threshold (dashed, blue)
    damaged_cell_umi <- 10^(em_params[["umi_means"]][which.max(em_params[["nf_means"]])])
    cell_umi <- 10^(em_params[["umi_means"]][which.min(em_params[["nf_means"]])])
    umi_diff <- log10(cell_umi - cell_umi*(umi_thresh/100))
    p3 <- p3 + ggplot2::geom_vline(xintercept = umi_diff, linetype="dashed", colour = cell_colour)
  }

  p.grid <- ggpubr::ggarrange(ggpubr::ggarrange(p1, p2, p3, ncol = 3))
  return(p.grid)

}
custom_identify_damaged_cells <- function(nf_umi_ed_ct,
                                   nf_sep=0.15,
                                   umi_sep_perc=50, # UMI counts percentage less than cell
                                   output_plots=FALSE,
                                   verbose=TRUE){

  # Check nf_umi_ed_ct argument
  if (any(class(nf_umi_ed_ct) == "data.frame")) {

    # Check four columns exist
    if(ncol(nf_umi_ed_ct)!=4){
      stop(paste0("nf_umi_ed_ct should be a data frame with four columns, see function arguments"), call.=FALSE)
    }

    # Assume nuclear fraction is in the first column
    nf <- unlist(nf_umi_ed_ct[, 1], use.names = FALSE)
    # Assume UMI counts are in the second column
    umi <- unlist(nf_umi_ed_ct[, 2], use.names = FALSE)
    # Assume third column contains "cell" or "empty_droplet"
    ed <- unlist(nf_umi_ed_ct[, 3], use.names = FALSE)
    # Assume fourth column contains cell types
    ct <- unlist(nf_umi_ed_ct[, 4], use.names = FALSE)


    # Check values are reasonable
    if(any(c(max(nf)>1, min(nf)<0))){
      warning(paste0("The nuclear fraction values provided in the first column of 'nf_umi_ed_ct' should be between 0 and 1, but values outside this range were identified : minimum = ",min(nf),", maximum = ",max(nf)), call.=FALSE)
    }

    if(!all(umi == floor(umi))){
      non_integer_examples <- which(umi != floor(umi))
      if(length(non_integer_examples)>5){
        non_integer_examples <- non_integer_examples[1:5]
      }
      non_integer_examples <- paste(umi[non_integer_examples], collapse = ",")
      warning(paste0("Non-integer values detected in the second column of 'nf_umi_ed_ct' (e.g. ",non_integer_examples,") where umi counts were expected"), call.=FALSE)
    }

    if(max(umi)<100){
      warning(paste0("The total umi counts provided in the second column of 'nf_umi_ed_ct' appear to be quite low (max = ",max(umi),"), are these the total UMI counts per cell?"), call.=FALSE)
    }

    if(!all(unique(ed)%in%c("cell", "empty_droplet"))){
      ed_output <- unique(ed)
      if(length(ed_output)>5){
        ed_output <- ed_output[1:5]
        ed_output <- paste(ed_output, collapse = ",")
      }
      warning(paste0("The third column of 'nf_umi_ed_ct' was expected to contain either 'cell' or 'empty_droplet' but contains; ",ed_output), call.=FALSE)
    }

    if(verbose){
      ct <- unique(ct)
      ct <- paste(ct, collapse = ",")
      print(paste0("The following cell types were provided; ", ct))
    }

  } else {
    stop(paste0("A data frame should be supplied to the nf_umi_ed_ct argument, but an object of class ",paste(class(nf_umi), collapse = "/")," was provided"), call.=FALSE)
  }

  # Extract data for EM
  em.data <- data.frame(nf = unlist(nf_umi_ed_ct[,1], use.names = FALSE),
                        umi = log10(unlist(nf_umi_ed_ct[,2], use.names = FALSE)),
                        ct = unlist(nf_umi_ed_ct[,4], use.names = FALSE))
  row.names(em.data) <- 1:nrow(em.data)

  # Filter out any empty droplets
  em.data <- em.data[nf_umi_ed_ct[,3]=="cell",]

  # Split by cell type
  em.data.ct <-split(em.data, em.data$ct)
  library(mclust)

  # Run EM for all cell types
  if(verbose){ print("Fitting models with EM") }
  em_mods <- lapply(em.data.ct, function(x) mclust::Mclust(data = x[,1:2], G = 1:2, modelNames = "EEI", verbose = verbose))

  # Assign barcodes as "cell" or "damaged_cell" using the `assess_EM` function
  em_mods_assessed  <- lapply(em_mods,
                              assess_EM,
                              nf_thresh = nf_sep,
                              umi_thresh = umi_sep_perc)

  # Create plots if requested
  if(output_plots){
    if(verbose){ print("Creating requested plots") }
    em_plots <- lapply(seq_along(em_mods), function(x) plot_EM(em = em_mods[[x]],
                                                               em_c = em_mods_assessed[[x]],
                                                               input_cell_type = names(em_mods)[x],
                                                               umi_thresh = umi_sep_perc,
                                                               nf_thresh = nf_sep))
    names(em_plots) <- names(em_mods)
  }

  # Update the input data frame "nf_umi_ed_ct" (which still contains empty
  # droplets) any with damaged_cell info
  names(em_mods_assessed) <- NULL
  em_mods_assessed <- do.call(rbind, em_mods_assessed)
  nf_umi_ed_ct$cell_status[as.integer(row.names(em_mods_assessed))] <- em_mods_assessed$classification

  # If plots were not requested, return results as a data frame
  if(output_plots){
    return(list(df=nf_umi_ed_ct, plots=em_plots))
  } else {
    return(list(df=nf_umi_ed_ct, plots=NULL))
  }
}


# ********* Consistent coloring for the cell in all three plots ********************
tmp_dc <-
custom_identify_damaged_cells(
    nf_umi_ed_ct = as_data_frame(cbind(tmp_ed, (
      big_sc$predicted.celltype.l2
    ))),
    nf_sep = 0.05,
    umi_sep_perc = 15,
    output_plots = TRUE
  )

# ************************************
library(ggplot2)
## Saving the determined "damaged cell" plots per annotated cell types in the dataset
lapply(names(tmp_dc[[2]]), function(x) {
  ggsave(paste0(x, "_nuc_frac_based_damaged_cells_plot.png"), tmp_dc[[2]][[x]], device = "png", dpi = 100,
         path = "docs/nuc_frac_plots/",height = 4,width = 12)
})

big_sc$nuc_frac <-metdat$nuc_frac
metdat$cell_status <- tmp_dc$df$cell_status
big_sc$cell_status <- tmp_dc$df$cell_status
qs::qsave(big_sc,file = "data/non_qc-trimmed_nf-based_cell_status_addedd_big_sc.qs")
