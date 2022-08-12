# Imports.
library(Seurat)
library(Signac)
library(MuDataSeurat)

################################### Liu ########################################

file_path <- "Phd/mowgli_reproducibility/data/Liu/liu_preprocessed.h5mu.gz"
out_path <- "Phd/mowgli_reproducibility/data/Liu/liu_seurat.RDS"

# Read the preprocessed MuData file as a Seurat object.
seurat_object <- MuDataSeurat::ReadH5MU(file_path)

#Scale and center the RNA data, then perform PCA.
Seurat::DefaultAssay(seurat_object) <- "rna"
VariableFeatures(seurat_object) <- rownames(seurat_object[["rna"]])
seurat_object <- Seurat::ScaleData(seurat_object)
seurat_object <- Seurat::RunPCA(seurat_object)

# Define all peaks features as top peaks, then run SVD.
Seurat::DefaultAssay(seurat_object) <- "atac"
seurat_object <- Signac::FindTopFeatures(seurat_object, min.cutoff = NULL)
seurat_object <- Signac::RunSVD(seurat_object)

# Run Seurat for Liu and simulated. We don't do tf-idf, so atac should not be called 'lsi'.
seurat_object <- Seurat::FindMultiModalNeighbors(
  seurat_object,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:50, 2:50),
  knn.range = 50
)

# Save Seurat object as a RDS.
saveRDS(seurat_object, out_path)

################################## Sim 1 #######################################

file_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_1.h5mu.gz"
out_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_1_seurat.RDS"

# Read the preprocessed MuData file as a Seurat object.
seurat_object <- MuDataSeurat::ReadH5MU(file_path)

#Scale and center the RNA data, then perform PCA.
Seurat::DefaultAssay(seurat_object) <- "rna"
VariableFeatures(seurat_object) <- rownames(seurat_object[["rna"]])
seurat_object <- Seurat::ScaleData(seurat_object)
seurat_object <- Seurat::RunPCA(seurat_object)

# Define all peaks features as top peaks, then run SVD.
Seurat::DefaultAssay(seurat_object) <- "atac"
seurat_object <- Signac::FindTopFeatures(seurat_object, min.cutoff = NULL)
seurat_object <- Signac::RunSVD(seurat_object)

# Run Seurat for Liu and simulated. We don't do tf-idf, so atac should not be called 'lsi'.
seurat_object <- Seurat::FindMultiModalNeighbors(
  seurat_object,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:50, 2:50),
  knn.range = 50
)

# Save Seurat object as a RDS.
saveRDS(seurat_object, out_path)

################################## Sim 2 #######################################

file_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_2.h5mu.gz"
out_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_2_seurat.RDS"

# Read the preprocessed MuData file as a Seurat object.
seurat_object <- MuDataSeurat::ReadH5MU(file_path)

#Scale and center the RNA data, then perform PCA.
Seurat::DefaultAssay(seurat_object) <- "rna"
VariableFeatures(seurat_object) <- rownames(seurat_object[["rna"]])
seurat_object <- Seurat::ScaleData(seurat_object)
seurat_object <- Seurat::RunPCA(seurat_object)

# Define all peaks features as top peaks, then run SVD.
Seurat::DefaultAssay(seurat_object) <- "atac"
seurat_object <- Signac::FindTopFeatures(seurat_object, min.cutoff = NULL)
seurat_object <- Signac::RunSVD(seurat_object)

# Run Seurat for Liu and simulated. We don't do tf-idf, so atac should not be called 'lsi'.
seurat_object <- Seurat::FindMultiModalNeighbors(
  seurat_object,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:50, 2:50),
  knn.range = 50
)

# Save Seurat object as a RDS.
saveRDS(seurat_object, out_path)

################################## Sim 3 #######################################

file_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_3.h5mu.gz"
out_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_3_seurat.RDS"

# Read the preprocessed MuData file as a Seurat object.
seurat_object <- MuDataSeurat::ReadH5MU(file_path)

#Scale and center the RNA data, then perform PCA.
Seurat::DefaultAssay(seurat_object) <- "rna"
VariableFeatures(seurat_object) <- rownames(seurat_object[["rna"]])
seurat_object <- Seurat::ScaleData(seurat_object)
seurat_object <- Seurat::RunPCA(seurat_object)

# Define all peaks features as top peaks, then run SVD.
Seurat::DefaultAssay(seurat_object) <- "atac"
seurat_object <- Signac::FindTopFeatures(seurat_object, min.cutoff = NULL)
seurat_object <- Signac::RunSVD(seurat_object)

# Run Seurat for Liu and simulated. We don't do tf-idf, so atac should not be called 'lsi'.
seurat_object <- Seurat::FindMultiModalNeighbors(
  seurat_object,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:50, 2:50),
  knn.range = 50
)

# Save Seurat object as a RDS.
saveRDS(seurat_object, out_path)

################################## Sim 4 #######################################

file_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_4.h5mu.gz"
out_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_4_seurat.RDS"

# Read the preprocessed MuData file as a Seurat object.
seurat_object <- MuDataSeurat::ReadH5MU(file_path)

#Scale and center the RNA data, then perform PCA.
Seurat::DefaultAssay(seurat_object) <- "rna"
VariableFeatures(seurat_object) <- rownames(seurat_object[["rna"]])
seurat_object <- Seurat::ScaleData(seurat_object)
seurat_object <- Seurat::RunPCA(seurat_object)

# Define all peaks features as top peaks, then run SVD.
Seurat::DefaultAssay(seurat_object) <- "atac"
seurat_object <- Signac::FindTopFeatures(seurat_object, min.cutoff = NULL)
seurat_object <- Signac::RunSVD(seurat_object)

# Run Seurat for Liu and simulated. We don't do tf-idf, so atac should not be called 'lsi'.
seurat_object <- Seurat::FindMultiModalNeighbors(
  seurat_object,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:50, 2:50),
  knn.range = 50
)

# Save Seurat object as a RDS.
saveRDS(seurat_object, out_path)

################################## Sim 5 #######################################

file_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_5.h5mu.gz"
out_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_5_seurat.RDS"

# Read the preprocessed MuData file as a Seurat object.
seurat_object <- MuDataSeurat::ReadH5MU(file_path)

#Scale and center the RNA data, then perform PCA.
Seurat::DefaultAssay(seurat_object) <- "rna"
VariableFeatures(seurat_object) <- rownames(seurat_object[["rna"]])
seurat_object <- Seurat::ScaleData(seurat_object)
seurat_object <- Seurat::RunPCA(seurat_object)

# Define all peaks features as top peaks, then run SVD.
Seurat::DefaultAssay(seurat_object) <- "atac"
seurat_object <- Signac::FindTopFeatures(seurat_object, min.cutoff = NULL)
seurat_object <- Signac::RunSVD(seurat_object)

# Run Seurat for Liu and simulated. We don't do tf-idf, so atac should not be called 'lsi'.
seurat_object <- Seurat::FindMultiModalNeighbors(
  seurat_object,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:50, 2:50),
  knn.range = 50
)

# Save Seurat object as a RDS.
saveRDS(seurat_object, out_path)

################################## PBMC ########################################

# Define the input and output paths.
file_path <- "Phd/mowgli_reproducibility/data/10X_PBMC_10k/pbmc_preprocessed.h5mu.gz"
out_path <- "Phd/mowgli_reproducibility/data/10X_PBMC_10k/pbmc_seurat.RDS"

# Read the preprocessed MuData file as a Seurat object.
seurat_object <- MuDataSeurat::ReadH5MU(file_path)

#Scale and center the RNA data, then perform PCA.
Seurat::DefaultAssay(seurat_object) <- "rna"
VariableFeatures(seurat_object) <- rownames(seurat_object[["rna"]])
seurat_object <- Seurat::ScaleData(seurat_object)
seurat_object <- Seurat::RunPCA(seurat_object)

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

# Save Seurat object as a RDS.
saveRDS(seurat_object, out_path)

################################## BMCITE ######################################

seurat_object <- MuDataSeurat::ReadH5AD("Phd/mowgli_reproducibility/data/BMCITE/bmcite_preprocessed_rna.h5ad")
seurat_object <- SeuratObject::RenameAssays(seurat_object, RNA = "rna")
seurat_object[["adt"]] <- MuDataSeurat::ReadH5AD("Phd/mowgli_reproducibility/data/BMCITE/bmcite_preprocessed_adt.h5ad")[['RNA']]
out_path <- "Phd/mowgli_reproducibility/data/BMCITE/bmcite_seurat.RDS"

#Scale and center the RNA data, then perform PCA.
Seurat::DefaultAssay(seurat_object) <- "rna"
VariableFeatures(seurat_object) <- rownames(seurat_object[["rna"]])
seurat_object <- Seurat::ScaleData(seurat_object)
seurat_object <- Seurat::RunPCA(seurat_object)

# Scale and center the ADT data, then perform PCA.
Seurat::DefaultAssay(seurat_object) <- "adt"
VariableFeatures(seurat_object) <- rownames(seurat_object[["adt"]])
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object, reduction.name = 'apca')

# Run Seurat for CITE-seq.
seurat_object <- Seurat::FindMultiModalNeighbors(
  seurat_object,
  reduction.list = list("pca", "apca"),
  dims.list = list(1:50, 1:24),
)

# Save Seurat object as a RDS.
saveRDS(seurat_object, out_path)

################################## OPCITE ######################################

seurat_object <- MuDataSeurat::ReadH5AD("Phd/mowgli_reproducibility/data/OPCITE/opcite_preprocessed_rna.h5ad")
seurat_object <- SeuratObject::RenameAssays(seurat_object, RNA = "rna")
seurat_object[["adt"]] <- MuDataSeurat::ReadH5AD("Phd/mowgli_reproducibility/data/OPCITE/opcite_preprocessed_adt.h5ad")[["RNA"]]
out_path <- "Phd/mowgli_reproducibility/data/OPCITE/opcite_seurat.RDS"

#Scale and center the RNA data, then perform PCA.
Seurat::DefaultAssay(seurat_object) <- "rna"
VariableFeatures(seurat_object) <- rownames(seurat_object[["rna"]])
seurat_object <- Seurat::ScaleData(seurat_object)
seurat_object <- Seurat::RunPCA(seurat_object)

# Scale and center the ADT data, then perform PCA.
Seurat::DefaultAssay(seurat_object) <- "adt"
VariableFeatures(seurat_object) <- rownames(seurat_object[["adt"]])
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object, reduction.name = 'apca')

# Run Seurat for CITE-seq.
seurat_object <- Seurat::FindMultiModalNeighbors(
  seurat_object,
  reduction.list = list("pca", "apca"),
  dims.list = list(1:50, 1:50),
)

# Save Seurat object as a RDS.
saveRDS(seurat_object, out_path)

################################ OPMULTIOME ####################################

seurat_object <- MuDataSeurat::ReadH5AD("Phd/mowgli_reproducibility/data/OP_multiome/opmultiome_preprocessed_rna.h5ad")
seurat_object <- SeuratObject::RenameAssays(seurat_object, RNA = "rna")
seurat_object[["atac"]] <- MuDataSeurat::ReadH5AD("Phd/mowgli_reproducibility/data/OP_multiome/opmultiome_preprocessed_atac.h5ad")[["RNA"]]
out_path <- "Phd/mowgli_reproducibility/data/OP_multiome/opmultiome_seurat.RDS"

#Scale and center the RNA data, then perform PCA.
Seurat::DefaultAssay(seurat_object) <- "rna"
VariableFeatures(seurat_object) <- rownames(seurat_object[["rna"]])
seurat_object <- Seurat::ScaleData(seurat_object)
seurat_object <- Seurat::RunPCA(seurat_object)

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

# Save Seurat object as a RDS.
saveRDS(seurat_object, out_path)

################################### TEA ########################################

file_path <- "Phd/mowgli_reproducibility/data/TEA/tea_preprocessed.h5mu.gz"
out_path <- "Phd/mowgli_reproducibility/data/TEA/tea_seurat.RDS"

# Read the preprocessed MuData file as a Seurat object.
seurat_object <- MuDataSeurat::ReadH5MU(file_path)

#Scale and center the RNA data, then perform PCA.
Seurat::DefaultAssay(seurat_object) <- "rna"
VariableFeatures(seurat_object) <- rownames(seurat_object[["rna"]])
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

# Run Seurat for TEA-seq. We don't do tf-idf, so atac should not be called 'lsi'.
seurat_object <- Seurat::FindMultiModalNeighbors(
  seurat_object,
  reduction.list = list("pca", "apca", "lsi"),
  dims.list = list(1:50, 1:45, 2:50),
)

# Save Seurat object as a RDS.
saveRDS(seurat_object, out_path)
