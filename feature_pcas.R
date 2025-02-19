library(ggplot2)
library(ggpubr)
library(ggsci)
library(umap)
library(plotly)
library(pandoc)

source("patient_values.R")


SPECIALS <- c("ENSG00000143546"="S100A8", "ENSG00000163220"="S100A9")


vlens <- function(x, n=2) {
  if (n == 2) {
    return(rowSums(abs(x)))
  }
  if (n == 2) {
    return(sqrt(rowSums(x^2)))
  }
  return((rowSums(x^n))^(1/n))
}


get_feature_type <- function(feature) {
  if (startsWith(feature, "Frequency")) return("scCluster")
  else if (grepl(feature, "__")) return("zOTU")
  return("BulkGene")
}


extract_timepoints <- function(data, metadata, x_col="PC1", y_col="PC2") {
  samples <- rownames(data)
  # column order: x_start, x_end, y_start, y_end
  timepoint_data <- data.frame(matrix(ncol = 3, nrow = 0))
  for (patient in unique(metadata[samples, "Patient"])) {
    patient_samples <- rownames(metadata[samples,])[
      metadata[samples, "Patient"] == patient]
    if (length(patient_samples) > 1) {
      pts <- metadata[patient_samples, "Timepoint"]
      names(pts) <- patient_samples
      pls <- metadata[patient_samples, "Location"]
      for (location in unique(pls)) {
        sub_samp <- rownames(metadata[patient_samples,])[pls == location]
        sub_samp <- sub_samp[order(pts[sub_samp])]
        # NOTE: performance-wise this is horrible (never use repeated rbind)
        #       but I'm too lazy to write it out the nice way and it works
        #       fine for the amount of data we have
        if (length(sub_samp) > 1) {
          timepoint_data <- rbind(
            timepoint_data,
            c(
              data[sub_samp[2], x_col], data[sub_samp[1], x_col],
              data[sub_samp[2], y_col], data[sub_samp[1], y_col]
            )
          )
          # check if "timepoint 3" exists
          if (length(sub_samp) == 3) {
            timepoint_data <- rbind(
              timepoint_data,
              c(
                data[sub_samp[3], x_col], data[sub_samp[2], x_col],
                data[sub_samp[3], y_col], data[sub_samp[2], y_col]
              )
            )
          }
        }
      }
    }
  }
  colnames(timepoint_data) <- c("x_start", "x_end", "y_start", "y_end")

  return(timepoint_data)
}


plot_embedding <- function(
  x, groups, ellipse=FALSE, locations=NULL, metadata=NULL, ...
) {
  if (!is.null(metadata)) {
    cols <- c("Patient", "Location", "NewNotInflamed", "Timepoint")
    x <- cbind(
      x, metadata[rownames(x), cols])
  }

  if (is.null(locations)) {
    data <- as.data.frame(x) %>%
      add_column(Sample=rownames(x)) %>%
      add_column(Group=groups[.$Sample])
    if (is.null(metadata)) {
      p <- ggplot(data, aes(x=PC1, y=PC2, colour=Group))
    } else {
      p <- ggplot(
        data,
        aes(
          x=PC1, y=PC2, colour=Group, text=paste0(
            "Patient: ", Patient, "<br>",
            "Location: ", Location, "<br>",
            "NotInflamedGroup: ", NewNotInflamed, "<br>",
            "Timepoint: ", Timepoint, "<br>"
          )
        )
      )
    }
  } else {
    data <- as.data.frame(x) %>%
      add_column(Sample=rownames(x)) %>%
      add_column(Group=groups[.$Sample])
    data$Location <- locations[data$Sample]
    if (is.null(metadata)) {
      p <- ggplot(data, aes(x=PC1, y=PC2, colour=Group, shape=Location))
    } else {
      p <- ggplot(
        data,
        aes(
          x=PC1, y=PC2, colour=Group, shape=Location, text=paste0(
            "Patient: ", Patient, "<br>",
            "Location: ", Location, "<br>",
            "NotInflamedGroup: ", NewNotInflamed, "<br>",
            "Timepoint: ", Timepoint, "<br>"
          )
        )
      )
    }
  }
  if (ellipse) {
    p <- p + stat_ellipse(aes(group=Group, fill=Group), alpha=.25, geom="polygon")
  }
  return(p + geom_point(...))
}


fit_pca_ <- function(data, features, samples=NULL) {
  if (is.null(samples)) {
    pca <- prcomp(data[,features], scale.=TRUE, center=TRUE)
  } else{
    sample_diff <- setdiff(samples, rownames(data))
    if (length(sample_diff) > 0) {
      warning(
        paste0("Some samples are not in the data: ", paste(sample_diff, ","))
      )
    }
    pca <- prcomp(
      data[intersect(samples, rownames(data)),features],
      scale.=TRUE, center=TRUE
    )
  }
  return(pca)
}


plot_pca <- function(
  data, groups, features, samples=NULL, ellipse=FALSE, locations=NULL,
  metadata=NULL, return_full=FALSE, ...
) {
  pca <- fit_pca_(data, features, samples)
  p <- plot_embedding(pca$x[,1:2], groups, ellipse, locations, metadata, ...)
  var_expl <- pca$sdev^2/sum(pca$sdev^2)
  p <- p + theme_pubr() +
    xlab(paste0('PC1 (', round(var_expl[1], 4)*100, '%)')) +
    ylab(paste0('PC2 (', round(var_expl[2], 4)*100, '%)'))

  if (return_full) {
    return(list(plot=p, pca=pca))
  }
  return(p)

}


plot_loadings <- function(
  data, features=NULL, samples=NULL, top_features=NULL, components=c(1, 2), ...
) {
  if (is.null(features)) features <- colnames(data)
  pca <- fit_pca_(data, features, samples)
  var_expl <- pca$sdev^2/sum(pca$sdev^2)

  rotations <- pca$rotation[,components]
  if (!is.null(top_features)) {
    lens <- vlens(rotations)
    top_feats <- rownames(pca$rotation)[order(lens, decreasing=TRUE)]
    rotations <- rotations[top_feats[1:top_features],]
  }
  rotations <- as.data.frame(rotations)
  rotations$DataType <- sapply(rownames(rotations), get_feature_type)

  ggplot(data=rotations, aes(x=0, y=0, xend=PC1, yend=PC2, colour=DataType)) +
    # arrow=arrow(length=unit(1, "picas"))
    geom_segment() +
    geom_text(
      aes(x=PC1, y=PC2, label=rownames(rotations), colour=DataType)) +
    theme_pubr() +
    xlab(paste0('PC1 (', round(var_expl[1], 4)*100, '%)')) +
    ylab(paste0('PC2 (', round(var_expl[2], 4)*100, '%)'))
}


plot_umap <- function(
  data, groups, features, samples=NULL, ellipse=FALSE, locations=NULL,
  metadata=NULL, ...
) {
  if (is.null(samples)) {
    embed <- umap(data[,features], scale.=TRUE, center=TRUE)
  } else{
    sample_diff <- setdiff(samples, rownames(data))
    if (length(sample_diff) > 0) {
      warning(
        paste0("Some samples are not in the data: ", paste(sample_diff, ","))
      )
    }
    embed <- umap(
      data[intersect(samples, rownames(data)),features],
      scale.=TRUE, center=TRUE
    )
  }

  # for compatibility with `plot_embedding`
  colnames(embed$layout) <- c("PC1", "PC2")
  p <- plot_embedding(embed$layout, groups, ellipse, locations, metadata, ...)
  return(
    p +  theme_pubr() +
      xlab("UMAP1") +
      ylab("UMAP2")
  )
}


to_interactive <- function(
  plot=NULL, plot_function=NULL, metadata=NULL, data=NULL, groups=NULL,
  features=NULL, samples=NULL, ellipse=FALSE, locations=NULL, ...
) {
  if (is.null(plot)) {
    plot <- plot_function(
      data, groups, features, samples, ellipse, locations, metadata, ...
    )
  }
  return(
    ggplotly(plot, tooltip=c("x", "y", "colour", "shape", "text"))
  )
}


cluster_labels <- function(data, features=NULL, samples=NULL, k=2) {
  if (is.null(features)) {
    features <- colnames(data)
  }
  if (is.null(samples)) {
    data <- data[,features]
  } else{
    sample_diff <- setdiff(samples, rownames(data))
    if (length(sample_diff) > 0) {
      warning(
        paste0("Some samples are not in the data: ", paste(sample_diff, ","))
      )
    }
    data <- data[intersect(samples, rownames(data)),features]
  }

  return(kmeans(data, k)$cluster)
}


# NOTE: this function only loads bulk-sc data
load_features <- function(
  file_base, target_group="both", file_ending="_no IBD_subnet.csv"
) {
  features <- c()
  if (target_group == "both") {
    for (group in c("Active", "Endoscopic_Healing")) {
      net_file <- paste0(file_base, group, file_ending)
      features <- unique(c(features, load_network_results(net_file)$nodes))
    }
  } else {
    net_file <- paste0(file_base, target_group, file_ending)
    features <- unique(load_network_results(net_file)$nodes)
  }

  # these are incorrectly named
  special_mask <- features %in% names(SPECIALS)
  if (any(special_mask)) {
    features[special_mask] <- SPECIALS[features[special_mask]]
  }

  return(features)
}


update_names <- function(annot, map, startstr, map_col="New_annotation") {
  mask <- startsWith(annot, startstr)
  annot[mask] <- map[annot[mask], map_col]
  return(annot)
}
