# Imports.
library(Seurat)
library(Signac)
library(MuDataSeurat)

# Define the input and output paths.
file_path <- "Phd/mowgli_reproducibility/data/10X_PBMC_10k/pbmc_preprocessed.h5mu.gz"
out_path <- "Phd/mowgli_reproducibility/data/10X_PBMC_10k/pbmc_seurat.RDS"

# Read the preprocessed MuData file as a Seurat object.
seurat_object <- MuDataSeurat::ReadH5MU(file_path)

# Scale and center the RNA data, then perform PCA.
Seurat::DefaultAssay(seurat_object) <- "rna"
seurat_object <- Seurat::ScaleData(seurat_object)
seurat_object <- Seurat::RunPCA(seurat_object)

# Define all peaks features as top peaks, then run PCA.
Seurat::DefaultAssay(seurat_object) <- "atac"
seurat_object <- Signac::FindTopFeatures(seurat_object, min.cutoff = NULL)
seurat_object <- Signac::RunSVD(seurat_object)

# Run Seurat. We don't do tf-idf, so atac should not be called 'lsi'.
seurat_object <- Seurat::FindMultiModalNeighbors(
  seurat_object,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:50, 2:50),
  knn.range = 200
)

# Save Seurat object as a RDS.
saveRDS(seurat_object, out_path)
