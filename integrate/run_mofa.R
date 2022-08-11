# Imports.
library(MOFA2)
library(MuDataSeurat)

# Define the hyperparameters.
num_factors <- 15
my_seed <- 42

# Define the input and output paths.
file_path <- "Phd/mowgli_reproducibility/data/10X_PBMC_10k/pbmc_preprocessed.h5mu.gz"
out_path <- paste0("Phd/mowgli_reproducibility/data/10X_PBMC_10k/pbmc_mofa_", num_factors, ".hdf5")

# Read the MuData file as a Seurat object.
seurat_object <- MuDataSeurat::ReadH5MU(file_path)

# Center the preprocessed RNA data.
Seurat::DefaultAssay(seurat_object) <- "rna"
seurat_object <- Seurat::ScaleData(seurat_object, do.center = TRUE, do.scale = FALSE)

# Create the mofa object from the seurat object.
mofa_object <- MOFA2::create_mofa(seurat_object, assays = c('rna', 'atac'))

# We don't need the Seurat object in memory anymore.
remove(seurat_object)

# Define the model options.
model_opts <- MOFA2::get_default_model_options(mofa_object)
model_opts$num_factors <- num_factors

# Define the training options.
train_opts <- MOFA2::get_default_training_options(mofa_object)
train_opts$seed <- my_seed

# Prepare the dimensionality reduction.
mofa_object <- MOFA2::prepare_mofa(
    object = mofa_object,
    model_options = model_opts,
    training_options = train_opts
)

# Perform the dimensionality reduction.
trained_mofa_object <- MOFA2::run_mofa(mofa_object, out_path, use_basilisk = FALSE)
