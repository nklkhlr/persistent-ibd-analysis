library(openxlsx)

DISTANCE_METRIC <- "euclidean"


extract_patient_distances <- function(
  dist_mat, patients, groups, locations, patient_subset
) {
  dense_dist <- as.matrix(dist_mat)

  # emtpy matrices for centroid distances
  centroid_distances <- list(
    Ileum=matrix(
      nrow=length(patient_subset), ncol=3,
      dimnames=list(patient_subset, c("noIBD_Active", "noIBD_EH", "Active_EH"))
    ),
    Sigma=matrix(
      nrow=length(patient_subset), ncol=3,
      dimnames=list(patient_subset, c("noIBD_Active", "noIBD_EH", "Active_EH"))
    )
  )
  # emtpy matrix for sample-wise distances
  sample_distances <- matrix(
    nrow=length(patient_subset), ncol=2,
    dimnames=list(patient_subset, c("Ileum", "Sigma"))
  )

  for (patient in patient_subset) {
    samples_i <- names(patients)[patients == patient]

    eh_samples <- samples_i[groups[samples_i] == "Endoscopic_Healing"]
    active_samples <- samples_i[groups[samples_i] == "Active"]

    patient <- as.character(patient)
    for (location in unique(locations[samples_i])) {
      location_eh <- eh_samples[locations[eh_samples] == location]
      location_active <- active_samples[locations[active_samples] == location]

      sample_distances[patient, location] <- mean(
        dense_dist[location_eh, location_active], na.rm=TRUE)

      # active to noibd centroid
      if (length(location_active) > 1) {
        centroid_distances[[location]][patient, "noIBD_Active"] <- mean(
          dense_dist[paste(location, "no IBD", sep="_"), location_active],
          na.rm=TRUE
        )
      } else if (length(location_active) == 1) {
        centroid_distances[[location]][patient, "noIBD_Active"] <- dense_dist[
          paste(location, "no IBD", sep="_"), location_active]
      }
      # EH to noibd and active centroid
      if (length(location_eh) > 1) {
        centroid_distances[[location]][patient, "noIBD_EH"] <- mean(
          dense_dist[paste(location, "no IBD", sep="_"), location_eh],
          na.rm=TRUE
        )
        centroid_distances[[location]][patient, "Active_EH"] <- rowMeans(
          dense_dist[paste(location, "Active", sep="_"), location_eh],
          na.rm=TRUE
        )
      } else if (length(location_eh) == 1){
        centroid_distances[[location]][patient, "noIBD_EH"] <- dense_dist[
          paste(location, "no IBD", sep="_"), location_eh]
        centroid_distances[[location]][patient, "Active_EH"] <- dense_dist[
          paste(location, "Active", sep="_"), location_eh]
      }
    }
  }
  return(
    list(
      centroids=centroid_distances,
      samples=sample_distances
    )
  )
}


plot_centroid_distances <- function(centroid_dists, ...) {
  add_column(as.data.frame(centroid_dists$Sigma), Location="Sigma") %>%
    rbind(add_column(as.data.frame(centroid_dists$Ileum), Location="Ileum")) %>%
    add_column(Patient=c(sapply(centroid_dists, rownames))) %>%
    pivot_longer(
       -c("Location", "Patient"), names_to="Type", values_to="Distance") %>%
    ggplot(aes(x=Type, y=Distance, label=Patient, color=Patient)) +
      geom_text(...) +
      facet_wrap(Location~Type, scales="free") +
      ylab(paste(DISTANCE_METRIC, "distance")) +
      scale_color_d3(palette="category20") +
      theme_pubr() +
      theme(axis.text.x=element_blank())
}


plot_sample_distances <- function(sample_dists, ...) {
  add_column(as.data.frame(sample_dists), Patient=rownames(sample_dists)) %>%
    pivot_longer(-"Patient", values_to="Distance", names_to="Location") %>%
    ggplot(aes(x=Location, y=Distance, label=Patient, color=Patient)) +
      geom_text(...) +
      facet_wrap(.~Location, scales="free") +
      ylab(paste(DISTANCE_METRIC, "distance")) +
      scale_color_d3(palette="category20") +
      theme_pubr() +
      theme(axis.text.x=element_blank())
}


save_distances <- function(distances, file) {
  sheet_list <- list(
    "SampleDistances"=distances$samples,
    "CentroidDistancesIleum"=distances$centroids$Ileum,
    "CentroidDistancesSigma"=distances$centroids$Sigma
  )
  write.xlsx(sheet_list, file, rowNames=TRUE)
}


cluster_from_distances <- function(distance_data, main_group=NULL) {
  if (is.list(distance_data)) {
    dists <- cbind(distance_data$Ileum, distance_data$Sigma)
    dists <- dists[,grepl(main_group, colnames(dists))]
  } else {
    dists <- distance_data
  }

  na_mask <- rowSums(is.na(dists)) == 0
  return(kmeans(dists[na_mask,], 2))
}

