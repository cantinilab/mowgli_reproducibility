import muon as mu
import torch
from torchnmf.nmf import NMF
import os
from rich.console import Console
import numpy as np
import sys

# Define the data and figure folder.
data_folder = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/"

# Check that we have a command line argument.
assert len(sys.argv) > 1

# Define data and figure paths for different datasets.
if sys.argv[1] == "bmcite":
    prefix = "BMCITE/bmcite"
    data_path = os.path.join(data_folder, prefix + "_preprocessed.h5mu.gz")

elif sys.argv[1] == "liu":
    prefix = "Liu/liu"
    data_path = os.path.join(data_folder, prefix + "_preprocessed.h5mu.gz")

elif sys.argv[1] == "sim1":
    prefix = "Liu/liu_simulated_1"
    data_path = os.path.join(data_folder, prefix + ".h5mu.gz")

elif sys.argv[1] == "sim2":
    prefix = "Liu/liu_simulated_2"
    data_path = os.path.join(data_folder, prefix + ".h5mu.gz")

elif sys.argv[1] == "sim3":
    prefix = "Liu/liu_simulated_3"
    data_path = os.path.join(data_folder, prefix + ".h5mu.gz")

elif sys.argv[1] == "sim4":
    prefix = "Liu/liu_simulated_4"
    data_path = os.path.join(data_folder, prefix + ".h5mu.gz")

elif sys.argv[1] == "sim5":
    prefix = "Liu/liu_simulated_5"
    data_path = os.path.join(data_folder, prefix + ".h5mu.gz")

elif sys.argv[1] == "opcite":
    prefix = "OPCITE/opcite"
    data_path = os.path.join(data_folder, prefix + "_preprocessed.h5mu.gz")

elif sys.argv[1] == "opmultiome":
    prefix = "OP_multiome/opmultiome"
    data_path = os.path.join(data_folder, prefix + "_preprocessed.h5mu.gz")

elif sys.argv[1] == "pbmc":
    prefix = "10X_PBMC_10k/pbmc"
    data_path = os.path.join(data_folder, prefix + "_preprocessed.h5mu.gz")

elif sys.argv[1] == "tea":
    prefix = "TEA/tea"
    data_path = os.path.join(data_folder, prefix + "_preprocessed.h5mu.gz")


console = Console()
with console.status("[bold green]Loading data...") as status:

    # Load the data.
    mdata = mu.read_h5mu(data_path)
    console.log("Data loaded.")

console = Console()
with console.status("[bold green]Performing NMF...") as status:

    # Define the rank.
    rank = 50

    # Define the NMF model.
    model = NMF(mdata.shape, rank=rank)

    # Fit the NMF model.
    model.fit(torch.hstack([torch.Tensor(mdata[mod].X) for mod in mdata.mod]), beta=2)

    # Write the embedding as an obsm.
    mdata.obsm["X_nmf"] = model.H.detach().numpy()

    # Save to disk.
    np.save(
        os.path.join(data_folder, prefix + f"_nmf_{rank}.npy"),
        mdata.obsm["X_nmf"],
    )

    console.log("NMF performed.")
