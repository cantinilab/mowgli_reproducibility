library(Seurat)
library(Signac)
library(MuDataSeurat)

file_path -> "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/10X_PBMC_10k/pbmc_preprocessed.h5mu.gz"
outfile -> "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/10X_PBMC_10k/pbmc_seurat.Rds"

seurat_object <- MuDataSeurat::ReadH5MU(file_path)

Seurat::DefaultAssay(seurat_object) <- "rna"

seurat_object <- Seurat::SCTransform(
    seurat_object,
    verbose = FALSE,
    assay = "rna"
)

seurat_object <- Seurat::RunPCA(seurat_object)

Seurat::DefaultAssay(seurat_object) <- "atac"

# Define all features as being top features.
seurat_object <- FindTopFeatures(seurat_object, min.cutoff = NULL)

# Run SVD.
seurat_object <- Signac::RunSVD(seurat_object)

# Run Seurat.
seurat_object <- Seurat::FindMultiModalNeighbors(
    seurat_object,
    reduction.list = list("pca", "lsi"),
    dims.list = list(1:50, 2:50),
    knn.range = 200
)

# Save seurat object.
writeRDS(seurat_object, outfile)