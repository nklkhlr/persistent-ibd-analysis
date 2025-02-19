source("revision_scoring.R")


#' Predict IBrD values for new samples
#'
#' Embeds new samples into the PCA space of the reference data.
#' Subsequently, samples are projected onto the principal curve
#' that was originally used to compute pseudotimes and IBrD values.
#'
#'
#' @param sce_obj `SingleCellExperiment` object containing the `slingshot`
#'                curve fitted on the original data
#' @param pca_obj `prcomp` object containing the PCA computed on the original 
#'                data. Must include centering/scaling parameters as well
#' @param bulk_data Normalized samples x gene matrix with bulk data of new 
#'                  samples to predict. Must contain all genes in `signature_genes`
#'                  as columns.
#' @param signature_genes Vector of genes in the IBrD signature
predict_new_samples <- function(
    sce_obj, 
    pca_obj, 
    bulk_data, 
    signature_genes
) {
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
    pca_data$new_scores <- abs(pca_data$new_times - max_time) / max_time

    return pca_data
}
