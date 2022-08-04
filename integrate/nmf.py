import muon as mu
import torch
from torchnmf.nmf import NMF
import os

# Define the path to the data.
data_path = os.path.join(
    "/users/csb/huizing/Documents/PhD/Code/",
    "mowgli_reproducibility/data/10X_PBMC_10k/",
    "pbmc_preprocessed.h5mu.gz"
)

# Load the data.
mdata = mu.read_h5mu(data_path)

# Define the NMF model.
model = NMF(mdata.shape, rank=30)

# Fit the NMF model.
model.fit(torch.hstack([torch.Tensor(mdata[mod].X) for mod in mdata.mod]), beta=2)

# Write the embedding as an obsm.
mdata.obsm["X_nmf"] = model.H.detach().numpy()