library(readxl)
library(tidyverse)
library(ggsci)

source("selected_scores.R")


# helper functions
extract_attr <- function(samples, sc_obj, col) {
    stopifnot(
        "'col' must be in 'sc_obj@colData'"=col %in% colnames(sc_obj@colData))
    
    attr <- sapply(
        samples,
        function(x) {
            tryCatch(
                sc_obj@colData[x, col],
                error=function(e) {return(NaN)}
            )
        }
    ) 
    names(attr) <- samples
    return(attr)
}


gene_annotation <- function(gene_data) {
    annots <- read_xlsx("Data/overlap.active.noIBD.EH.noIBD.xlsx") %>%
        column_to_rownames("ensID")
    gene_data <- gene_data[rownames(annots),]
    rownames(gene_data) <- annots[["label"]]
    
    return(gene_data)
}


if (sys.nframe() == 0) {
    load("revision/sce_pca.RData")

    # load and format data
    metadata <- read_excel("Data/Metadata_Activityscore.xlsx") %>%
        add_column(
            Score=extract_attr(.$`#SampleID_16S`, sce, "activity_score"),
            Pseudotime=extract_attr(.$`#SampleID_16S`, sce, "slingPseudotime_1")
        ) %>%
        as.data.frame
    rownames(metadata) <- metadata$Sample_ID_Transkriptome

    sig_genes <- read.csv("correlation_plots/features.csv", row.names=1)[,1]
    rna_data <- read_excel("Data/vst.batch.corrected.xlsx") %>%
        column_to_rownames("gene") %>%
        gene_annotation %>%
        t
    rownames(rna_data) <- metadata[rownames(rna_data), "#SampleID_16S"]
    rownames(metadata) <- metadata$`#SampleID_16S`

    # load PCA computed in `selected_scores` ('active_sig_plot$pca')
    # and compute PC coordinates of 'new' samples and plot together with
    pca_data <- to_pc_coords(rna_data[,sig_genes], active_sig_plot$pca) %>%
        as.data.frame %>%
        bind_cols(metadata[rownames(rna_data),])

    # compute embedded curve to plot
    embedded <- embedCurves(sce, "PCA")
    embedded <- slingCurves(embedded)[[1]]
    embedded <- data.frame(embedded$s[embedded$ord,])
     

    plot_unscored <- ggplot(
            pca_data,
            aes(
                x=PC1, y=PC2, color=Score, shape=SCORE_Notinflamed,
            )
        ) +
        geom_point(size=3) +
        geom_path(data=embedded, mapping=aes(x=PC1, y=PC2), inherit.aes=FALSE) +
        scale_color_viridis_c()
    ggsave(
        plot=plot_unscored, file="revision/unscored.pdf",
        width=16, height=9, device="pdf"
    )

    # Project samples onto curve to get pseudotimes/scores for new samples
    n_pcs <- get_n_pcs(active_sig_plot$pca)
    curve <- slingCurves(sce)$Lineage1

    new_projection <- princurve::project_to_curve(
        as.matrix(pca_data[,sapply(1:n_pcs, function(i) paste0("PC", i))]),
        s=curve$s[curve$ord, , drop=FALSE],
        stretch=0
    )
    new_projection$lambda <- new_projection$lambda - min(new_projection$lambda)
    pca_data$new_times <- new_projection$lambda[rownames(pca_data)]

    # pseudotime to activity score
    max_time <- max(slingPseudotime(sce))
    pca_data$new_scores <- abs(pca_data$new_times - max_time) / max_time

    pca_data$is_old <- is.na(pca_data$Pseudotime)

    # plotting
    rescored_new <- ggplot(
            pca_data,
            aes(x=PC1, y=PC2, shape=is_old, color=new_scores)
        ) +
        geom_point(size=3) +
        geom_path(data=embedded, mapping=aes(x=PC1, y=PC2), inherit.aes=FALSE) +
        scale_color_viridis_c() +
        theme_minimal()
    ggsave(
        plot=rescored_new, file="revision/rescored_set_annotation.pdf",
        width=16, height=9, device="pdf"
    )

    rescored_status <- ggplot(
            pca_data,
            aes(x=PC1, y=PC2, shape=SCORE_Notinflamed, color=new_scores)
        ) +
        geom_point(size=3) +
        geom_path(data=embedded, mapping=aes(x=PC1, y=PC2), inherit.aes=FALSE) +
        scale_color_viridis_c() +
        theme_minimal()
    ggsave(
        plot=rescored_status, file="revision/rescored_disease_status_large.pdf",
        width=16, height=9, device="pdf"
    )

    write.csv(
        data.frame(new=pca_data$new_times, old=pca_data$Pseudotime),
        file="revision/inferred_pseudotimes.csv"
    )

    write.csv(
        data.frame(new=pca_data$new_scores, old=sce@colData[rownames(pca_data), "activity_score"]),
        file="revision/inferred_scores.csv"
    )
}
