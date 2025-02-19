source("feature_pcas.R")
source("distances.R")

library(tidyverse)
library(ggbeeswarm)
library(see)
library(SingleCellExperiment)
library(slingshot)
library(glmnet)
library(readxl)


load_features <- function(file) {
  corrs <- read.csv(file=file, row.names=1, check.names=FALSE)
  return(c(rownames(corrs), colnames(corrs)))
}

load_bulk_data <- function(file_path, id_map) {
  bulk_data <- read_xlsx(file_path) %>%
    column_to_rownames("gene")
  colnames(bulk_data) <- id_map[
    gsub("^X", "", colnames(bulk_data))]
  bulk_data <- bulk_data[,!is.na(colnames(bulk_data))]
  return(bulk_data)
}

col_to_vec <- function(df, col) {
  vec <- df[,col]
  names(vec) <- rownames(df)
  return(vec)
}

set_dist_names <- function(dmat, sample_names) {
  mat <- as.matrix(dmat)
  dimnames(mat) <- list(sample_names, sample_names)
  return(as.dist(mat))
}


compute_centroids <- function(
  sce_obj, groups, locations, col="activity_score"
) {
  patient_mask <- names(groups) %in% rownames(sce_obj@colData)
  ileum_group_centroids <- sapply(
    c("no IBD", "Active"),
    function(group) {
      group_mask <- groups == group & locations == "Ileum"
      group_samples <- names(groups)[patient_mask & group_mask]
      mean(sce_obj@colData[group_samples, col], na.rm=TRUE)
    }
  )
  names(ileum_group_centroids) <- paste0("Ileum_", names(ileum_group_centroids))

  sigma_group_centroids <- sapply(
    c("no IBD", "Active"),
    function(group) {
      group_mask <- groups == group & locations == "Sigma"
      group_samples <- names(groups)[patient_mask & group_mask]
      mean(sce_obj@colData[group_samples, col], na.rm=TRUE)
    }
  )
  names(sigma_group_centroids) <- paste0("Sigma_", names(sigma_group_centroids))

  return(
    list(
      Ileum=ileum_group_centroids,
      Sigma=sigma_group_centroids
    )
  )
}


compute_distances <- function(
  sce_obj, patient_ids, groups, locations, eh_patients, col="activity_score"
) {
  centroids <- compute_centroids(sce_obj, groups, locations, col)

  snames <- c(rownames(sce_obj@colData), c(sapply(centroids, names)))
  dist_data <- c(sce_obj@colData[,col], centroids$Ileum, centroids$Sigma)
  dists <- dist(dist_data, "manhattan") %>%
    set_dist_names(snames)

  return(
    extract_patient_distances(
      dists, patient_ids, groups, locations, eh_patients)
  )
}


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


prepare_plot_data_ <- function(exp_obj, groups, shape_name, color_by_activity) {
  plot_data <- as.data.frame(reducedDim(exp_obj, "PCA"))
  plot_data[shape_name] <- groups[rownames(plot_data)]
  name_order <- rownames(plot_data)
  if (color_by_activity) {
    plot_data$ActivityScore <- exp_obj@colData[name_order, "activity_score"]
  } else {
    plot_data$Pseudotime <- exp_obj@colData[name_order, "slingPseudotime_1"]
  }
  return(plot_data)
}


plot_pca_pseudotime <- function(
  exp_obj, emb_curve, groups, color_by_activity=TRUE, shape_name="NotInflamed",
  interactive_annotations=NULL, ...
) {
  if (color_by_activity) {
    color <- "ActivityScore"
  } else {
    color <- "Pseudotime"
  }
  plot_data <- prepare_plot_data_(
    exp_obj, groups, shape_name, color_by_activity)

  if (is.null(interactive_annotations)) {
    p <- ggplot(
      plot_data, aes(x=PC1, y=PC2, color=!!sym(color), shape=!!sym(shape_name))
    )
  } else {
    p <- add_column(
      plot_data,
      Patient=interactive_annotations[rownames(plot_data), "Patient"],
      Location=interactive_annotations[rownames(plot_data), "Location"],
      Disease=interactive_annotations[rownames(plot_data), "Disease"]
    ) %>%
      ggplot(
        aes(
          x=PC1, y=PC2, color=!!sym(color), shape=!!sym(shape_name),
          text=paste0(
            "Patient: ", Patient, "<br>",
            "Location: ", Location, "<br>",
            "Disease: ", Disease, "<br>"
          )
        ),
      )
  }

  return(
    p +
      geom_point(...) +
      scale_color_viridis_c() +
      geom_path(
        data=emb_curve, aes(x=PC1, y=PC2), linewidth=1.2, inherit.aes=FALSE)
  )
}


invert_times_ <- function(
  col_data, groups, control_group="no IBD", active_group="Active"
) {
  min_ <- min(col_data$slingPseudotime_1)
  min_sample <- rownames(col_data)[which(col_data$slingPseudotime_1 == min_)]
  if (is.vector(min_sample)) return(all(groups[min_sample] != control_group))
  return(groups[min_sample] == active_group)
}


activity_score <- function(sce_object, groups, ...) {
  col_frame <- sce_object@colData
  if (invert_times_(col_frame, groups, ...)) {
    pseudotime <- abs(
      col_frame$slingPseudotime_1 - max(col_frame$slingPseudotime_1))
  } else {
    pseudotime <- col_frame$slingPseudotime_1
  }

  max_time <- max(pseudotime)
  sce_object$activity_score <- pseudotime / max_time
  return(sce_object)
}


plot_score_distributions <- function(
  exp_obj, groups, swarmplot=TRUE, subsamples=NULL
) {
  activity <- exp_obj@colData$activity_score
  names(activity) <- rownames(exp_obj@colData)

  if (!is.null(subsamples)) {
    activity <- activity[subsamples]
  }

  p <- data.frame(ActivityScore=activity, Group=groups[names(activity)]) %>%
    ggplot(aes(x=Group, y=ActivityScore, color=Group, fill=Group))
  if (swarmplot) return(p + geom_beeswarm(cex=3))
  return(p + geom_violindot())
}


get_eh_origin <- function(metadata) {
  eh_samples <- rownames(metadata)[metadata$NewNotInflamed == "Endoscopic_Healing"]
  patients <- col_to_vec(metadata[eh_samples,], "Patient")
  locations <- col_to_vec(metadata[eh_samples,], "Location")
  timepoints <- col_to_vec(metadata[eh_samples,], "Timepoint")
  mode(timepoints) <- "integer"

  predecessors <- sapply(
    eh_samples,
    function(sample) {
      if (timepoints[sample] == 1) {
        # if EH is first timepoint we cannot find a predecessor
        warning(
          paste0(
            "'", sample,
            "' is the first timepoint of this patient. No predecessor available!"
          )
        )
        return(NA)
      }

      c("Location", "Timepoint", "NewNotInflamed")

      patient_mask <- metadata$Patient == patients[sample]
      location_mask <- metadata$Location == locations[sample]
      if (sum(patient_mask & location_mask) == 1) {
        # this means we can't find the same patient and location in any state
        # other than EH
        warning(
          paste0(
            "'", sample,
            "' is the only sample of this patient from the ", locations[sample],
            ". No predecessor available!"
          )
        )
        return(NA)
      }

      sub_time <- col_to_vec(
        metadata[patient_mask & location_mask,], "Timepoint")
      mode(sub_time) <- "integer"

      predecessor <- names(sub_time)[sub_time == sub_time[sample] - 1]
      if (is.null(predecessor)) return(NA)

      return(metadata[predecessor, "NewNotInflamed"])
    }
  )
  return(predecessors)
}


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


#' Predict pseudotime for new data points
#'
#' @param new_data N x 69 matrix with new data points to compute pseudotime for
#'   or Principal Components computed with the original PCA object
#' @param sling_obj `SingleCellExperiment` containing `slingshot` object to
#'  predict pseudotimes
#' @param pca PCA object with which original data was embedded. Must be given
#'   if `new_data` contains the original data points and must not be given if
#'   `new_data` contains the principal components
#'
#' @return a `PseudotimeOrdering` object containing pseudotimes for new data
predict_pseudotime <- function(new_data, sce_sling_obj, groups, pca=NULL) {
  if (!is.null(pca)) {
    if (ncol(new_data) != length(pca$center)) {
      stop(
        paste0(
          "`new_data` must contain exactly ", length(pca$center), " features (",
          "i.e. columns). If `new_data` contains new samples embedded in the ",
          "original PC space make sure to leave `pca` NULL"
        )
      )
    }
    new_pca <- to_pc_coords(new_data, pca)
  }
  else
    new_pca <- new_data

  browser()
  n_pcs <- ncol(reducedDim(sce_sling_obj, "PCA"))
  ptime <- slingshot::predict(sce_sling_obj$slingshot, new_pca[,1:n_pcs])

  new_times <- data.frame(pathStat(ptime))
  colnames(new_times) <- c("slingPseudotime_1")
  max_time <- max(sce_sling_obj@colData$slingPseudotime_1)
  if (invert_times_(sce_sling_obj@colData, groups)) {
    new_times$activity_score <- abs(new_times$slingPseudotime_1 - max_time) / max_time
  } else{
    new_times$activity_score <- new_times$slingPseudotime_1 / max_time
  }

  new_sce <- SingleCellExperiment(
    t(new_pca), reducedDims=SimpleList(PCA=new_pca[,1:n_pcs]))
  new_sce@colData <- S4Vectors::DataFrame(new_times)

  return(new_sce)
}


summary_plots <- function(
  sce_obj, emb_curve, groups, annotations, file_base, 
  show_interactive=TRUE, store_with_metadata=TRUE, ...
) {
  pca_plot <- plot_pca_pseudotime(sce_obj, emb_curve, groups, ...) +
    theme_bw()
  ggsave(
    filename=paste0(file_base, "PCA.pdf"), plot=pca_plot, width=16, height=9,
    device="pdf"
  )

  if (show_interactive) {
    plot_pca_pseudotime(
      sce_obj, emb_curve, groups, interactive_annotations=annotations, ...
    ) %>%
      plotly::ggplotly() %>%
      htmlwidgets::saveWidget(paste0(file_base, "PCA.html"))
  }
  
  swarm_plot <- plot_score_distributions(sce_obj, groups) +
    theme_bw()
  ggsave(
    filename=paste0(file_base, "swarm.pdf"), plot=swarm_plot, width=16,
    height=9, device="pdf"
  )

  violin_plot <- plot_score_distributions(sce_obj, groups, FALSE) +
    theme_bw()
  ggsave(
    filename=paste0(file_base, "violin.pdf"), plot=violin_plot, width=16,
    height=9, device="pdf"
  )

  if (store_with_metadata) {
    add_column(
      as.data.frame(sce_obj@colData[,2:3]),
      Disease=disease[rownames(sce_obj@colData)],
      Location=locations[rownames(sce_obj@colData)],
      NewOverall=metadata[rownames(sce_obj@colData), "NewOverall"],
      NewNotInflamed=not_inflamed[rownames(sce_obj@colData)]
    ) %>%
      write.csv(file=paste0(file_base, "Scores.csv"), row.names=TRUE)
  } else {
    write.csv(
      as.data.frame(sce_obj@colData[,2:3]),
      file=paste0(file_base, "Scores.csv"), row.names=TRUE
    )
  }
}


#' Compute effect sizes of Treatment and Individuals over all features
#' individually while correcting for Patient-specific effects with a linear 
#' mixed model of the form 
#' `Feature ~ (intercept) + Group + Treatment + (1|Patient)`
#'
#' @param num_data Feature data with samples in rows. All columns will be tested
#'                 individually
#' @param meta_data `data.frame` containing Diet, Activity and Patient annotation
#' @param intercept Whether to fit the LMM with an intercept
#' @param ... Optional parameters to pass to `nlme::lmer`
#'
#' @returns list of length m with per-feature model outputs
pca_lmm <- function(
  num_data, meta_data, intercept=TRUE, with_disease=FALSE, ...
) {
  samples <- intersect(rownames(num_data), rownames(meta_data)) 
  lm_data <- cbind(
    num_data[samples,],
    meta_data[samples, c("NewNotInflamed", "Location", "Patient", "AntiTNF", "Disease")]
  )
  
  if (intercept) 
    base_formula <- " ~ 1 + Location + AntiTNF + (1|Patient)"
  else 
    base_formula <- " ~ Location + AntiTNF + (1|Patient)" 
  
  if (with_disease)
    base_formula <- paste0(base_formula, " + Disease")
  else
    base_formula <- paste0(base_formula, " + NewNotInflamed")
  
  models <- lapply(
    colnames(num_data),
    function(feature) {
      lm_out <- lmerTest::lmer(
        as.formula(paste0(feature, base_formula)), data=lm_data, ...)
    }
  ) 
  names(models) <- colnames(num_data)
  
  return(models)
}


if (sys.nframe() == 0) {
  library(argparse)

  parser <- ArgumentParser()
  parser$add_argument(
    "--path", type="character", default=paste0(getwd(), "/correlation_plots"),
    help="Folder where correlation network files are located"
  )
  parser$add_argument(
    "--exclude_sc", type="logical", default=TRUE,
    help="Whether to exclude single-cell data; Defaults to TRUE"
  )
  parser$add_argument(
    "--remove-not-inflamed", type="logical", default=TRUE,
    help="Whether to exclude 'not inflamed' samples; Defaults to TRUE"
  )
  args <- parser$parse_args()
  path <- args$path

  load("Data/sc_processed.RData")
  load("Data/metadata.RData")

  # ======== #
  ### main ###
  # ======== #
  # making sure required dirs exist
  dir.create(paste0(path, "/ActivityScore"), recursive=TRUE)
  dir.create(paste0(path, "/plots"), recursive=TRUE)

  # metadata setup
  new_metadata <- readxl::read_excel("Data/marker_metadata") %>%
    column_to_rownames("#SampleID_16S")
  metadata <- add_column(
    metadata,
    NewOverall=new_metadata[rownames(metadata), "SCORE_Overall"],
    NewNotInflamed=new_metadata[rownames(metadata), "SCORE_Notinflamed"],
    BulkID=new_metadata[rownames(metadata), "Sample_ID_Transkriptome"]
  ) %>%
    dplyr::rename(Patient=CED)
  names(metadata)[names(metadata) == "Anti TNF"] <- "AntiTNF"

  # extracting certain metadata
  not_inflamed <- metadata$NewNotInflamed
  names(not_inflamed) <- rownames(metadata)

  locations <- metadata$Location
  names(locations) <- rownames(metadata)

  disease <- metadata$Disease
  names(disease) <- rownames(metadata)


  # load bulk data
  bulk_id_map <- rownames(metadata)
  names(bulk_id_map) <- metadata$BulkID

  bulk_data <- load_bulk_data("Data/vst.batch.corrected.xlsx", bulk_id_map)

  bulk_annot <- read_xlsx("Data/overlap.active.noIBD.EH.noIBD.xlsx") %>%
    column_to_rownames("ensID")

  bulk_data <- bulk_data[rownames(bulk_annot),]
  rownames(bulk_data) <- bulk_annot[["label"]]

  # combine into one df
  combined_samples <- colnames(bulk_data) %>%
    intersect(rownames(sc_cluster))

  combined_data <- cbind.data.frame(
    sc_cluster[combined_samples,], sc_subcluster[combined_samples,],
    t(bulk_data[,combined_samples])
  )

  # load significant features
  if (args$exclude_sc) {
    active_sig <- read.csv(paste0(path, "/selected_genes.csv"), header=FALSE) %>%
      t %>%
      as.vector
  }
  else{ 
    active_sig <- load_features(paste0(path, "/sc_bulk_network.csv"))
    active_sig <- active_sig[!grepl("^Frequency", active_sig)]
  }
  write.csv(active_sig, file=paste0(path, "/features.csv"))

  # remove not_inflamed samples
  if (args$remove_not_inflamed) {
    sub_samples <- names(not_inflamed)[not_inflamed != "Not_inflamed"]
  } else{
    sub_samples <- names(not_inflamed)
  }

  # perform pca
  active_sig_plot <- plot_pca(
    combined_data, not_inflamed, active_sig, samples=sub_samples,
    ellipse=TRUE, metadata=metadata, return_full=TRUE, size=3
  )
  arrow_plot <- active_sig_plot$plot +
    geom_segment(
      aes(x=x_end, y=y_end, xend=x_start, yend=y_start),
      data=extract_timepoints(active_sig_plot$pca$x, metadata),
      inherit.aes=FALSE, arrow=arrow(length=unit(0.01, "npc"))
    )
  ggsave(
    filename=paste0(path, "/plots/PCA_Active_features_arrows.pdf"), plot=arrow_plot,
    height=9, width=16, device="pdf"
  )
  to_interactive(arrow_plot) %>%
    as_widget %>%
    htmlwidgets::saveWidget(
      paste0(path, "/plots/PCA_Active_features_arrows.html"), selfcontained=TRUE)

  active_sig_loadings <- plot_loadings(
    combined_data, active_sig, samples=sub_samples)
  ggsave(
    filename=paste0(path, "/plots/PCA_Active_features_loadings.pdf"),
    plot=active_sig_loadings,
    height=9, width=16, device="pdf"
  )
  
  save(
    active_sig_plot, locations, disease, not_inflamed,
    file=paste0(path, "/trajectory_data.RData")
  )
  
  # lmm testing for principal components
  metadata$NewNotInflamed <- factor(
    metadata$NewNotInflamed,
    levels=c("no IBD", "Active", "Endoscopic_Healing", "Not_inflamed")
  )
  pca_lmms <- pca_lmm(active_sig_plot$pca$x, metadata)
  pvals <- sapply(pca_lmms, function(x) summary(x)$coefficients[,"Pr(>|t|)"])
  qvals <- t(apply(pvals, 1, p.adjust, method="bonferroni"))
  write.csv(qvals[,1:10], file=paste0(path, "/pca_testing.csv"))
  
  disease_samples <- rownames(active_sig_plot$pca$x)[
    metadata[rownames(active_sig_plot$pca$x), "Disease"] != "no IBD"]
  pca_lmms <- pca_lmm(active_sig_plot$pca$x, metadata, with_disease=TRUE)
  pvals <- sapply(pca_lmms, function(x) summary(x)$coefficients[,"Pr(>|t|)"])
  qvals <- t(apply(pvals, 1, p.adjust, method="bonferroni"))
  write.csv(qvals[,1:10], file=paste0(path, "/pca_testing_non_ibd_no_disease.csv"))

  # ============== #
  ### score main ###
  # ============== #
  n_pcs <- get_n_pcs(active_sig_plot$pca)
  cat("Number of PCs explaining 80% variance: ", n_pcs, "\n")
  sce <- SingleCellExperiment(
    t(active_sig_plot$pca$x),
    reducedDims=SimpleList(PCA=active_sig_plot$pca$x[,1:n_pcs])
  )

  sce <- slingshot(sce, reducedDim="PCA")
  sce <- activity_score(sce, not_inflamed)

  embedded <- embedCurves(sce, "PCA")
  embedded <- slingCurves(embedded)[[1]]
  embedded <- data.frame(embedded$s[embedded$ord,])
 
  # ======================================== # 
  ##### Testing of new sample prediction #####
  # ======================================== #
  new_old_samples <- combined_data[1:4, active_sig]
  predict_pseudotime(new_old_samples, sce, not_inflamed[1:4], active_sig_plot$pca) 

  # plotting results
  file_base <- paste0(path, "/ActivityScore/ActivityScore_")
  summary_plots(sce, embedded, not_inflamed, metadata, file_base, size=5)

  location_plot <- plot_pca_pseudotime(
    sce, embedded, locations, shape_name="Location", size=5) +
    stat_ellipse() +
    theme_bw()
  ggsave(
    filename=paste0(path, "/ActivityScore/ActivityScore_PCA_Locations.pdf"),
    plot=location_plot, width=16, height=9, device="pdf"
  )
  disease_plot <- plot_pca_pseudotime(
    sce, embedded, disease, shape_name="Disease", size=5) +
    stat_ellipse() +
    theme_bw()
  ggsave(
    filename=paste0(path, "/ActivityScore/ActivityScore_PCA_Disease.pdf"),
    plot=disease_plot, width=16, height=9, device="pdf"
  )

  # ===================== #
  # score-based distances #
  # ===================== #
  patients <- read_excel(paste0(path, "/PatientPlots/Active_Distances.xlsx"))
  patients <- c(patients, metadata$Patient[metadata$NewOverall == "no IBD"])
  patient_samples <- rownames(metadata)[metadata$Patient %in% patients] %>%
    intersect(rownames(active_sig_plot$pca$x))

  groups <- metadata[patient_samples, "NewOverall"]
  names(groups) <- patient_samples
  locations <- metadata[patient_samples, "Location"]
  names(locations) <- patient_samples
  patient_ids <- metadata[patient_samples, "Patient"]
  names(patient_ids) <- patient_samples

  eh_annotations <- read_excel("Data/EH_Group.xlsx") %>% rename(Patient=CED)
  eh_samples <- rownames(metadata)[metadata$Patient %in% eh_annotations$Patient]
  eh_patients <- unique(metadata[eh_samples,"Patient"])

  eh_patient_dists <- compute_distances(
    sce, patient_ids, groups, locations, eh_patients)
  save_distances(eh_patient_dists, paste0(path, "/ActivityScore/Distances.xlsx"))

  plot_centroid_distances(eh_patient_dists$centroids, size=6) %>%
    ggsave(
      filename=paste0(path, "/ActivityScore/Distances_Centroid_EH.pdf"), plot=.,
      width=16, height=16, device="pdf"
    )
  plot_sample_distances(eh_patient_dists$samples, size=6) %>%
    ggsave(
      filename=paste0(path, "/ActivityScore/Distances_Sample_EH.pdf"), plot=.,
      width=16, height=16, device="pdf"
    )

  eh_origin <- get_eh_origin(metadata)
  eh_origin_nona <- eh_origin[!is.na(eh_origin)]
  plot_score_distributions(
    sce, eh_origin_nona, subsamples=names(eh_origin_nona)) %>%
    ggsave(
      filename=paste0(path, "/ActivityScore/PredecessorScores.pdf"), plot=.,
      width=16, height=9, device="pdf"
    )
}
