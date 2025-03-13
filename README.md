# Code to for the IBrD computations in [Siebert et al.]()


## Prerequisites

The code in this repository was tested with R version 4.2.1 
For required `R` packages and the used versions please see `sessionInfo.txt`

## Adding new samples

To compute the IBrD value for new samples, use the function `predict_new_samples` in `predict_new_samples.R`.

```R
# loads two objects named "od_pca" (prcomp object) and "sce" (SingleCellExperiment object)
load("sce_pca.RData")

# load new (bulk RNAseq) data
new_data <- read.csv("data/new_data.csv", row.names=1)

# predict IBrD value for new samples
updated_pca <- predict_new_samples(od_pca$pca, sce, new_data)
# plot values for old and new samples as colors and indicate new samples
# with different shape
plot_pca_pseudotime(sce, updated_pca, "new_IBrD", "is_old")
```

To obtain the original `SCE` and `prcomp` objects, please contact the lead authors.

