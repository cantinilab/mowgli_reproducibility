######################################## IMPORTS #########################################

# Load libraries.
import scanpy as sc
import muon as mu
from rich.console import Console
import os
import numpy as np
import sys

########################################## PATHS #########################################

# Define the data and figure folder.
data_folder = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/"
figure_folder = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/img/"

# Check that we have a command line argument.
assert(len(sys.argv) > 1)

# Define data and figure paths for different datasets.
if sys.argv[1] == "bmcite":
    data_path = os.path.join(data_folder, "BMCITE/bmcite_preprocessed.h5mu.gz")
    figure_path = os.path.join(figure_folder, "bmcite/")

elif sys.argv[1] == "liu":
    data_path = os.path.join(data_folder, "Liu/liu_preprocessed.h5mu.gz")
    figure_path = os.path.join(figure_folder, "liu/")

elif sys.argv[1] == "sim1":
    data_path = os.path.join(data_folder, "Liu/liu_simulated_1.h5mu.gz")
    figure_path = os.path.join(figure_folder, "sim1/")

elif sys.argv[1] == "sim2":
    data_path = os.path.join(data_folder, "Liu/liu_simulated_2.h5mu.gz")
    figure_path = os.path.join(figure_folder, "sim2/")

elif sys.argv[1] == "sim3":
    data_path = os.path.join(data_folder, "Liu/liu_simulated_3.h5mu.gz")
    figure_path = os.path.join(figure_folder, "sim3/")

elif sys.argv[1] == "sim4":
    data_path = os.path.join(data_folder, "Liu/liu_simulated_4.h5mu.gz")
    figure_path = os.path.join(figure_folder, "sim4/")

elif sys.argv[1] == "sim5":
    data_path = os.path.join(data_folder, "Liu/liu_simulated_5.h5mu.gz")
    figure_path = os.path.join(figure_folder, "sim5/")

elif sys.argv[1] == "opcite":
    data_path = os.path.join(data_folder, "OPCITE/opcite_preprocessed.h5mu.gz")
    figure_path = os.path.join(figure_folder, "opcite/")

elif sys.argv[1] == "opmultiome":
    data_path = os.path.join(data_folder, "OP_multiome/opmultiome_preprocessed.h5mu.gz")
    figure_path = os.path.join(figure_folder, "opmultiome/")

elif sys.argv[1] == "pbmc":
    data_path = os.path.join(data_folder, "10X_PBMC_10k/pbmc_preprocessed.h5mu.gz")
    figure_path = os.path.join(figure_folder, "pbmc/")

elif sys.argv[1] == "tea":
    data_path = os.path.join(data_folder, "TEA/tea_preprocessed.h5mu.gz")
    figure_path = os.path.join(figure_folder, "tea/")


##################################### LOAD DATA ##########################################

console = Console()
with console.status("[bold green]Loading data...") as status:

    # Load the original data.
    mdata = mu.read_h5mu(data_path)
    console.log("Data loaded.")

#################################### PLOT UMAPS ##########################################

with console.status("[bold green]Plotting UMAPs...") as status:

    # Define figure path.
    sc.settings.figdir = figure_path

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
