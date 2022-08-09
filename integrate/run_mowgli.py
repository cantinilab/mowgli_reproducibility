from mowgli import models
import numpy as np
import muon as mu
import os

# Define the path to the data.
data_path = os.path.join(
    "/users/csb/huizing/Documents/PhD/Code/",
    "mowgli_reproducibility/data/10X_PBMC_10k/",
    "pbmc_preprocessed.h5mu.gz",
)

# Load the data.
mdata = mu.read_h5mu(data_path)

# Define the Mowgli model.
# model = models.MowgliModel(latent_dim=30, rho_h=1e-2, rho_w=1e-3)

# Fit the Mowgli model.
# model.train(mdata)

mdata.obsm["W_OT"] = np.load(
    os.path.join(
        "/users/csb/huizing/Documents/PhD/Code/Mowgli/local_analysis/from_jz/w/",
        "pbmc_mowgli_correlation_50_0_05_0_01_0_001.npy",
    ),
    allow_pickle=True,
).item()["W"]

# Save to disk.
np.save(
    os.path.join(
        "/users/csb/huizing/Documents/PhD/Code/",
        "mowgli_reproducibility/data/10X_PBMC_10k/pbmc_mowgli_50.npy",
    ),
    mdata.obsm["W_OT"],
)
