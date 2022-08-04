library(SeuratDisk)
library(MOFA2)
library(MuDataSeurat)

file_path -> "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/10X_PBMC_10k/pbmc_preprocessed.h5mu.gz"
outfile -> "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/10X_PBMC_10k/pbmc_mofa.hdf5"

# Create Seurat object from h5mu file.
seurat_object <- ReadH5MU(file_path)

# Create MOFA object from Seurat object.
mofa_object <- create_mofa_from_matrix(
    data = list(
        seurat_object@assays$rna@counts,
        seurat_object@assays$atac@counts
    )
)

# Define hyperparameters.
num_factors <- 30
my_seed <- 42

# Define data options.
data_opts <- get_default_data_options(mofa_object)

# Define model options.
model_opts <- get_default_model_options(mofa_object)
model_opts$num_factors <- num_factors

# Define training options.
train_opts <- get_default_training_options(mofa_object)
train_opts$seed <- my_seed

# Perform dimensionality reduction.
mofa_object <- prepare_mofa(
    object = mofa_object,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
)

# Save the MOFA object.
trained_mofa_object <- run_mofa(mofa_object, out_file)
