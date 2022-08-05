######################################## IMPORTS #########################################

# Load libraries.
import scanpy as sc
import muon as mu
from rich.console import Console
import os
import numpy as np

##################################### LOAD DATA ##########################################

console = Console()
with console.status("[bold green]Loading data...") as status:

    # Define the data path.
    data_path = os.path.join(
        "/users/csb/huizing/Documents/PhD/Code/",
        "mowgli_reproducibility/data/10X_PBMC_10k/",
        "pbmc_preprocessed.h5mu.gz",
    )

    # Load the original data.
    mdata = mu.read_h5mu(data_path)
    console.log("Data loaded.")

#################################### PLOT UMAPS ##########################################

with console.status("[bold green]Plotting UMAPs...") as status:

    # Define figure path.
    sc.settings.figdir = os.path.join(
        "/users/csb/huizing/Documents/PhD/Code/",
        "mowgli_reproducibility/img/pbmc/",
    )

    # Iterate over modalities.
    for mod in mdata.mod:

        # Compute the kNN graph.
        sc.pp.neighbors(mdata[mod])

        # Compute the UMAP embedding.
        sc.tl.umap(mdata[mod])

        # Define the figure name.
        figure_name = f"_{mod}.pdf"

        # Plot the UMAP.
        sc.pl.umap(
            mdata[mod],
            color="celltype",
            save=figure_name,
            show=False,
        )

        # Print the figure name.
        console.log(f"{figure_name} plotted.")

    # Load the NMF embedding.
    nmf_path = os.path.join(
        "/users/csb/huizing/Documents/PhD/Code/",
        "mowgli_reproducibility/data/10X_PBMC_10k/pbmc_nmf.npy",
    )
    mdata.obsm["X_nmf"] = np.load(nmf_path)
    mdata.uns = {}

    # Compute the kNN graph.
    sc.pp.neighbors(mdata, use_rep="X_nmf")

    # Compute the UMAP embedding.
    sc.tl.umap(mdata)

    # Define the figure name.
    figure_name = f"_nmf.pdf"

    # Plot the UMAP.
    sc.pl.umap(
        mdata,
        color="rna:celltype",
        save=figure_name,
        show=False,
    )

    # Print the figure name.
    console.log(f"{figure_name} plotted.")
