# Imports.
library(Seurat)
library(Signac)
library(MuData)

remotes::install_github("grimbough/rhdf5")
remotes::install_github("ilia-kats/MuData")

# Define the input and output paths.
file_path <- "Phd/mowgli_reproducibility/data/10X_PBMC_10k/pbmc_preprocessed.h5mu.gz"
out_path <- "Phd/mowgli_reproducibility/data/10X_PBMC_10k/pbmc_seurat.RDS"

file_path <- "Phd/mowgli_reproducibility/data/BMCITE/bmcite_preprocessed.h5mu.gz"
out_path <- "Phd/mowgli_reproducibility/data/BMCITE/bmcite_seurat.RDS"

file_path <- "Phd/mowgli_reproducibility/data/Liu/liu_preprocessed.h5mu.gz"
out_path <- "Phd/mowgli_reproducibility/data/Liu/liu_seurat.RDS"

file_path <- "Phd/mowgli_reproducibility/data/OP_multiome/opmultiome_preprocessed.h5mu.gz"
out_path <- "Phd/mowgli_reproducibility/data/OP_multiome/opmultiome_seurat.RDS"

file_path <- "Phd/mowgli_reproducibility/data/OPCITE/opcite_preprocessed.h5mu.gz"
out_path <- "Phd/mowgli_reproducibility/data/OPCITE/opcite_seurat.RDS"

file_path <- "Phd/mowgli_reproducibility/data/TEA/tea_preprocessed.h5mu.gz"
out_path <- "Phd/mowgli_reproducibility/data/TEA/tea_seurat.RDS"

file_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_1.h5mu.gz"
out_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_1_seurat.RDS"

file_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_2.h5mu.gz"
out_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_2_seurat.RDS"

file_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_3.h5mu.gz"
out_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_3_seurat.RDS"

file_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_4.h5mu.gz"
out_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_4_seurat.RDS"

file_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_5.h5mu.gz"
out_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_5_seurat.RDS"

# Read the preprocessed MuData file as a Seurat object.
seurat_object <- MuDataSeurat::ReadH5MU(file_path)

#Scale and center the RNA data, then perform PCA.
Seurat::DefaultAssay(seurat_object) <- "rna"
seurat_object <- Seurat::ScaleData(seurat_object)
seurat_object <- Seurat::RunPCA(seurat_object)

# Scale and center the ADT data, then perform PCA.
Seurat::DefaultAssay(seurat_object) <- "adt"
VariableFeatures(seurat_object) <- rownames(seurat_object[["adt"]])
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object, reduction.name = 'apca')

# Define all peaks features as top peaks, then run SVD.
Seurat::DefaultAssay(seurat_object) <- "atac"
seurat_object <- Signac::FindTopFeatures(seurat_object, min.cutoff = NULL)
seurat_object <- Signac::RunSVD(seurat_object)

# Run Seurat for RNA and ATAC. We don't do tf-idf, so atac should not be called 'lsi'.
seurat_object <- Seurat::FindMultiModalNeighbors(
  seurat_object,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:50, 2:50),
)

# Run Seurat for TEA-seq. We don't do tf-idf, so atac should not be called 'lsi'.
seurat_object <- Seurat::FindMultiModalNeighbors(
  seurat_object,
  reduction.list = list("pca", "apca", "lsi"),
  dims.list = list(1:50, 1:45, 2:50),
)

# Run Seurat for CITE-seq.
seurat_object <- Seurat::FindMultiModalNeighbors(
  seurat_object,
  reduction.list = list("pca", "apca"),
  dims.list = list(1:50, 1:50),
)

# Run Seurat for Liu and simulated. We don't do tf-idf, so atac should not be called 'lsi'.
seurat_object <- Seurat::FindMultiModalNeighbors(
  seurat_object,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:50, 2:50),
  knn.range = 50
)

# Save Seurat object as a RDS.
saveRDS(seurat_object, out_path)

