library(tidyverse)
library(readxl)
library(pheatmap)


# ========= #
# Functions #
# ========= #
load_network_results <- function(file) {
  edges <- read.csv(file=file, row.names=1)
  nodes <- unique(c(edges[,1], edges[,2]))
  return(list(nodes=nodes, edges=edges))
}


load_maaslin_results <- function(file, param, ctrl_group, group=NULL) {
  lmm_data <- read_tsv(file) %>%
    filter(metadata == param) %>%
    dplyr::select(c("feature", "coef", "value"))

  if (is.null(group)) {
    lmm_data$Association <- sapply(
      seq_len(nrow(lmm_data)),
      function(i) {
        if (lmm_data$coef[i] > 0) return(lmm_data$value[i])
        return(ctrl_group)
      }
    )
    lmm_data <- dplyr::select(lmm_data, -value)
  } else {
    lmm_data <- lmm_data %>%
      filter(value == group) %>%
      add_column(Association=ifelse(.$coef > 0, group, ctrl_group))
  }

  return(lmm_data)
}


load_feature_names <- function() {
  zotu_annot <- read.csv("Data/taxonomy_extracted.csv")
  bulk_annot <- read_excel("Data/BulkEnsembleIDS_GenNames.xlsx")

  feature_names <- c(paste0("g__", zotu_annot$Genus), bulk_annot$New_annotation)
  names(feature_names) <- c(zotu_annot$X, bulk_annot$ensID)
  return(feature_names)
}


map_feature_name <-  function(feature, annots) {
  if (startsWith(feature, "Zotu")) {
    return(paste0(annots[feature], "(", feature, ")"))
  }
  if (!is.na(annots[feature]) & !is.null(annots[feature])) {
    return(annots[feature])
  }
  return(feature)
}


get_modality <- function(feature_names) {
  sapply(
    feature_names,
    function(x) {
      if (startsWith(x, "Zotu")) return("zOTU")
      else if (startsWith(x, "ENSG")) return("bulkRNA")
      return("scCluster")
    }
  )
}


format_plot_names <- function(feature, annots) {
  if (startsWith(feature, "Zotu")) {
    return(paste0(substr(annots[feature], start=1, stop=25), " (", feature, ")"))
  }
  if (!is.null(annots[feature])) {
    return(substr(annots[feature], start=1, stop=25))
  }
  return(feature)
}


patient_heatmap <- function(
  sample_data, metadata, net_file, lmm_files, name_map, param, ctrl_group,
  group=NULL, ...
) {
  network <- load_network_results(net_file)
  network$nodes <- sapply(network$nodes, map_feature_name, annots=name_map)
  names(network)$nodes <- NULL

  lmm_res <- bind_rows(
    lapply(
      lmm_files,
      function(file) {
        if (grepl("zotu", file)) {
          return(
            load_maaslin_results(
              file=file, param="NewOverall", group=group, ctrl_group=ctrl_group)
          )
        }
        return(
          load_maaslin_results(
            file=file, param=param, group=group, ctrl_group=ctrl_group))
      }
    ),
    .id="Modality"
  )
  rownames(lmm_res) <- make.unique(sapply(
    lmm_res$feature, map_feature_name, annots=name_map))

  if (is.null(group)) {
    sample_meta <- metadata
  } else {
    sample_meta <- filter(metadata, metadata[[param]] %in% c(group, ctrl_group))
  }
  samples <- intersect(rownames(sample_meta), rownames(sample_data))

  # heatmap with patients sorted by: param -> Patient -> Location
  # features sorted by association
  sample_order <- order(
    sample_meta[samples, param], sample_meta[samples, "Patient"],
    sample_meta[samples, "Location"]
  )
  feature_order <- order(
    lmm_res[network$nodes,] %>% .$Association,
    lmm_res[network$nodes,] %>% .$Modality
  )

  annot_row <- data.frame(
    Patient=sample_meta[samples[sample_order], "Patient"],
    Location=sample_meta[samples[sample_order], "Location"]
  )
  annot_row[[param]] <- sample_meta[samples[sample_order], param]
  rownames(annot_row) <- samples[sample_order]

  annot_col <- data.frame(
    Association=lmm_res[network$nodes[feature_order],] %>% .$Association,
    Modality=lmm_res[network$nodes[feature_order],] %>% .$Modality
  )
  rownames(annot_col) <- network$nodes[feature_order]

  pheatmap(
    sample_data[samples[sample_order], network$nodes[feature_order]],
    cluster_rows=FALSE, cluster_cols=FALSE, scale="column",
    annotation_row=annot_row, annotation_col=annot_col,
    ...
  )
}
