library(ggsci)
library(readxl)
library(slingshot)
library(tidyverse)


# helper function to get the minimum number of PCs explaining 
# 80% total variance
get_n_pcs <- function(pca_obj, min_var_expl=.8, min_pcs=2) {
  var_expl <- pca_obj$sdev^2/sum(pca_obj$sdev^2)
  total_var <- sum(var_expl[1:min_pcs])
  n <- min_pcs
  while(total_var < min_var_expl) {
    n <- n + 1
    total_var <- total_var + var_expl[n]
  }

  return(n)
}


# helper function to rotate new sample points into PCA space
# based on the original data
to_pc_coords <- function(new_data, pca) {
  scaled_data <- sapply(
    1:ncol(new_data),
	function(i) {
		(new_data[,i] - pca$center[i]) / pca$scale[i]
	}
  )
  new_pca <- as.matrix(scaled_data) %*% pca$rotation
  dimnames(new_pca) <- list(rownames(new_data), colnames(pca$rotation))

  return(new_pca)
}


#' Predict IBrD values for new samples
#'
#' Embeds new samples into the PCA space of the reference data.
#' Subsequently, samples are projected onto the principal curve
#' that was originally used to compute pseudotimes and IBrD values.
#'
#'
#' @param sce_obj `SingleCellExperiment` object containing the `slingshot`
#'                curve fitted on the original data. 
#' @param pca_obj `prcomp` object containing the PCA computed on the original 
#'                data. Must include centering/scaling parameters as well
#' @param bulk_data Normalized samples x gene matrix with bulk data of new 
#'                  samples to predict. Must contain all genes in `signature_genes`
#'                  as columns.
predict_new_samples <- function(
    sce_obj, 
    pca_obj, 
    bulk_data 
) {
    signature_genes <- rownames(pca_obj$pca$rotation)
    # project new samples into PCA space
    # `to_pc_coords` also centers/scales the data
    pca_data <- to_pc_coords(bulk_data[,signature_genes], pca_obj) %>%
        as.data.frame
    # get the number of PCs to use (minimum number of PCs with a cumulative
    # explained variance >= 80%)
    n_pcs <- get_n_pcs(pca)

    # extract original curve and project new samples onto it
    curve <- slingCurves(sce)$Lineage1
    new_projection <- princurve::project_to_curve(
        as.matrix(pca_data[,sapply(1:n_pcs, function(i) paste0("PC", i))]),
        s=curve$s[curve$ord, , drop=FALSE],
        stretch=0
    )

    # lambda to pseudotime
    new_projection$lambda <- new_projection$lambda - min(new_projection$lambda)
    pca_data$new_times <- new_projection$lambda[rownames(pca_data)]
    # pseudotime to IBrD score
    max_time <- max(slingPseudotime(sce))
    pca_data$new_IBrD <- abs(pca_data$new_times - max_time) / max_time
    pca_data$is_old <- is.na(pca_data$Pseudotime)

    return pca_data
}


#' Plot IBrD values in a PCA plot
#'
#' @param sce_obj `SingleCellExperiment` object containing the `slingshot`
#'                curve fitted on the original data.
#' @param pca_data Data frame containing PCA coordinates and values to 
#'                plot as colors/shapes
#' @param color_name Name of column in `pca_data` to plot as colors
#' @param shape_name Name of column in `pca_data` to plot as shapes (optional)
#' @param ... Optional arguments to pass to `geom_point`
#'
plot_pca_pseudotime <- function(
    sce_obj,
    pca_data,
    color_name,
    shape_name=NULL,
    ...
) {
    if (is.null(shape_name)) {
        p <- ggplot(
            pca_data,
            aes(x=PC1, y=PC2, color=color_name)
        )
    } else {
        p <- ggplot(
            pca_data,
            aes(x=PC1, y=PC2, color=color_name, shape=shape_name)
        )
    }

    emb_curve <- embedCurve(sce, "PCA")
    p +
        geom_point(...) +
        geom_path(data=emb_curve, mapping=aes(x=PC1, y=PC2), inherit.aes=FALSE) +
        scale_color_viridis_c() +
        theme_minimal()
}

