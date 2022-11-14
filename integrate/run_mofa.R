# Imports.
library(MOFA2)
library(MuDataSeurat)

# Define the hyperparameters.
num_factors <- 5
my_seed <- 42

################################### Liu ########################################

file_path <- "Phd/mowgli_reproducibility/data/Liu/liu_preprocessed.h5mu.gz"
out_path <- paste0("Phd/mowgli_reproducibility/data/Liu/liu_mofa_", num_factors, ".hdf5")

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

################################## Sim 1 #######################################

file_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_1.h5mu.gz"
out_path <- paste0("Phd/mowgli_reproducibility/data/Liu/liu_simulated_1_mofa_", num_factors, ".hdf5")

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

################################## Sim 2 #######################################

file_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_2.h5mu.gz"
out_path <- paste0("Phd/mowgli_reproducibility/data/Liu/liu_simulated_2_mofa_", num_factors, ".hdf5")

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

################################## Sim 3 #######################################

file_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_3.h5mu.gz"
out_path <- paste0("Phd/mowgli_reproducibility/data/Liu/liu_simulated_3_mofa_", num_factors, ".hdf5")

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

################################## Sim 4 #######################################

file_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_4.h5mu.gz"
out_path <- paste0("Phd/mowgli_reproducibility/data/Liu/liu_simulated_4_mofa_", num_factors, ".hdf5")

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

################################## Sim 5 #######################################

file_path <- "Phd/mowgli_reproducibility/data/Liu/liu_simulated_5.h5mu.gz"
out_path <- paste0("Phd/mowgli_reproducibility/data/Liu/liu_simulated_5_mofa_", num_factors, ".hdf5")

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

################################## PBMC ########################################

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

################################## BMCITE ######################################

seurat_object <- MuDataSeurat::ReadH5AD("Phd/mowgli_reproducibility/data/BMCITE/bmcite_preprocessed_rna.h5ad")
seurat_object <- SeuratObject::RenameAssays(seurat_object, RNA = "rna")
seurat_object[["adt"]] <- MuDataSeurat::ReadH5AD("Phd/mowgli_reproducibility/data/BMCITE/bmcite_preprocessed_adt.h5ad")[['RNA']]
out_path <- paste0("Phd/mowgli_reproducibility/data/BMCITE/bmcite_mofa_", num_factors, ".hdf5")

# Center the preprocessed RNA data.
Seurat::DefaultAssay(seurat_object) <- "rna"
seurat_object <- Seurat::ScaleData(seurat_object, do.center = TRUE, do.scale = FALSE)

Seurat::DefaultAssay(seurat_object) <- "rna"
VariableFeatures(seurat_object) <- rownames(seurat_object[["rna"]])

Seurat::DefaultAssay(seurat_object) <- "adt"
VariableFeatures(seurat_object) <- rownames(seurat_object[["adt"]])

# Create the mofa object from the seurat object.
mofa_object <- MOFA2::create_mofa(seurat_object, assays = c('rna', 'adt'))

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

################################## OPCITE ######################################

seurat_object <- MuDataSeurat::ReadH5AD("Phd/mowgli_reproducibility/data/OPCITE/opcite_preprocessed_rna.h5ad")
seurat_object <- SeuratObject::RenameAssays(seurat_object, RNA = "rna")
seurat_object[["adt"]] <- MuDataSeurat::ReadH5AD("Phd/mowgli_reproducibility/data/OPCITE/opcite_preprocessed_adt.h5ad")[["RNA"]]
out_path <- paste0("Phd/mowgli_reproducibility/data/OPCITE/opcite_mofa_", num_factors, ".hdf5")

# Center the preprocessed RNA data.
Seurat::DefaultAssay(seurat_object) <- "rna"
seurat_object <- Seurat::ScaleData(seurat_object, do.center = TRUE, do.scale = FALSE)

Seurat::DefaultAssay(seurat_object) <- "rna"
VariableFeatures(seurat_object) <- rownames(seurat_object[["rna"]])

Seurat::DefaultAssay(seurat_object) <- "adt"
VariableFeatures(seurat_object) <- rownames(seurat_object[["adt"]])

# Create the mofa object from the seurat object.
mofa_object <- MOFA2::create_mofa(seurat_object, assays = c('rna', 'adt'))

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

################################ OPMULTIOME ####################################

seurat_object <- MuDataSeurat::ReadH5AD("Phd/mowgli_reproducibility/data/OP_multiome/opmultiome_preprocessed_rna.h5ad")
seurat_object <- SeuratObject::RenameAssays(seurat_object, RNA = "rna")
seurat_object[["atac"]] <- MuDataSeurat::ReadH5AD("Phd/mowgli_reproducibility/data/OP_multiome/opmultiome_preprocessed_atac.h5ad")[["RNA"]]
out_path <- paste0("Phd/mowgli_reproducibility/data/OP_multiome/opmultiome_mofa_", num_factors, ".hdf5")

# Center the preprocessed RNA data.
Seurat::DefaultAssay(seurat_object) <- "rna"
seurat_object <- Seurat::ScaleData(seurat_object, do.center = TRUE, do.scale = FALSE)

Seurat::DefaultAssay(seurat_object) <- "rna"
VariableFeatures(seurat_object) <- rownames(seurat_object[["rna"]])

Seurat::DefaultAssay(seurat_object) <- "atac"
VariableFeatures(seurat_object) <- rownames(seurat_object[["atac"]])

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

################################### TEA ########################################

file_path <- "Phd/mowgli_reproducibility/data/TEA/tea_preprocessed.h5mu.gz"
out_path <- paste0("Phd/mowgli_reproducibility/data/TEA/tea_mofa_", num_factors, ".hdf5")

# Read the MuData file as a Seurat object.
seurat_object <- MuDataSeurat::ReadH5MU(file_path)

# Center the preprocessed RNA data.
Seurat::DefaultAssay(seurat_object) <- "rna"
seurat_object <- Seurat::ScaleData(seurat_object, do.center = TRUE, do.scale = FALSE)

# Create the mofa object from the seurat object.
mofa_object <- MOFA2::create_mofa(seurat_object, assays = c('rna', 'atac', 'adt'))

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

