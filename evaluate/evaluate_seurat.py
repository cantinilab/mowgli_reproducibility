######################################## IMPORTS #########################################

# Load libraries.
import muon as mu
from rich.console import Console
import os
import numpy as np
from mowgli import score
import pickle
import pandas as pd

console = Console()

# Define the data and figure folder.
data_folder = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/"

# Define data paths for different datasets.
data_path = {
    "bmcite_seurat": data_folder + "BMCITE/bmcite_preprocessed.h5mu.gz",
    "liu_seurat": data_folder + "Liu/liu_preprocessed.h5mu.gz",
    "sim1_seurat": data_folder + "Liu/liu_simulated_1.h5mu.gz",
    "sim2_seurat": data_folder + "Liu/liu_simulated_2.h5mu.gz",
    "sim3_seurat": data_folder + "Liu/liu_simulated_3.h5mu.gz",
    "sim4_seurat": data_folder + "Liu/liu_simulated_4.h5mu.gz",
    "sim5_seurat": data_folder + "Liu/liu_simulated_5.h5mu.gz",
    "opcite_seurat": data_folder + "OPCITE/opcite_preprocessed.h5mu.gz",
    "opmultiome_seurat": data_folder + "OP_multiome/opmultiome_preprocessed.h5mu.gz",
    "pbmc_seurat": data_folder + "10X_PBMC_10k/pbmc_preprocessed.h5mu.gz",
    # "tea_seurat": data_folder + "TEA/tea_preprocessed.h5mu.gz",
}

seurat_path = {
    "bmcite_seurat": data_folder + "BMCITE/bmcite_seurat",
    "liu_seurat": data_folder + "Liu/liu_seurat",
    "sim1_seurat": data_folder + "Liu/liu_simulated_1_seurat",
    "sim2_seurat": data_folder + "Liu/liu_simulated_2_seurat",
    "sim3_seurat": data_folder + "Liu/liu_simulated_3_seurat",
    "sim4_seurat": data_folder + "Liu/liu_simulated_4_seurat",
    "sim5_seurat": data_folder + "Liu/liu_simulated_5_seurat",
    "opcite_seurat": data_folder + "OPCITE/opcite_seurat",
    "opmultiome_seurat": data_folder + "OP_multiome/opmultiome_seurat",
    "pbmc_seurat": data_folder + "10X_PBMC_10k/pbmc_seurat",
    # "tea_seurat": data_folder + "TEA/tea_seurat",
}

###################################### JACCARD THING #####################################


# Define the Jaccard index.
def jaccard(a, b):
    inter = len(np.intersect1d(a, b))
    return inter / (len(a) + len(b) - inter)


# Function reweighting Scanpy's kNN based on the Jaccard index.
def jaccard_denoising(adata):

    # Iterate over all cells.
    for i in range(adata.n_obs):

        # Get nearest neighbors of cell i.
        _, idx_i = adata.obsp["distances"][i].nonzero()

        # Iterate over the nearest neighbors.
        for j in idx_i:

            # Get the nearest neighbors of cell j.
            _, idx_j = adata.obsp["distances"][j].nonzero()

            # Compute the Jaccard index.
            d = jaccard(idx_i, idx_j)

            # Reweight the kNN.
            adata.obsp["connectivities"][i, j] = d
            adata.obsp["connectivities"][j, i] = d


################################# EVALUATING SEURAT ########################################

with console.status("[bold green]Evaluating Seurat...") as status:

    # Intialize a dictionary for the scores.
    scores_dict = {}

    # Set the range of nearest neighbors.
    k_range = list(range(1, 30))

    # Set the range of resulotions.
    res_range = list(np.arange(0.1, 2.1, 0.1))

    previous_path = ""

    for xp_name in data_path:

        # Load the data.
        console.log(f"Loading data for {xp_name} [bold green]")
        if previous_path != data_path[xp_name]:
            mdata = mu.read_h5mu(data_path[xp_name])
            previous_path = data_path[xp_name]
        console.log("Data loaded.")

        # Initialise scores for this experiment.
        scores_dict[xp_name] = {}

        # Log the experiment name.
        console.log(f"Starting to compute scores for {xp_name} [bold green]")

        # Load the seurat knn.
        knn = pd.read_csv(seurat_path[xp_name] + "_knn.csv", index_col=0).to_numpy() - 1
        mdata.uns = {}

        console.log("Loaded the knn [bold green]")

        # Compute the purity score for varying k nearest neighbors.
        purity_scores = []
        for k in k_range:

            # Log the value of k.
            console.log(f"Computing purity score for k={k} [bold green]")

            # Compute the purity score.
            s = score.knn_purity_score(knn[:, :k], mdata.obs["rna:celltype"])
            purity_scores.append(s)

        scores_dict[xp_name]["Purity scores"] = purity_scores
        scores_dict[xp_name]["k range"] = k_range

        console.log("Computed the purity scores. Phew! [bold green]")

        # Load the clustering.
        leiden = pd.read_csv(seurat_path[xp_name] + "_clustering.csv", index_col=0).to_numpy()

        # Compute the ARI for varying resolutions of Leiden clustering.
        aris = []
        for i, res in enumerate(res_range):

            # Log the value of resolution.
            console.log(f"Computing ARI for resolution={res} [bold green]")

            # Compute the ARI.
            aris.append(score.ARI(mdata.obs["rna:celltype"], leiden[:, i]))

        scores_dict[xp_name]["ARIs"] = aris
        scores_dict[xp_name]["res_range"] = res_range

        console.log("Computed the ARIs. Phew! [bold green]")

    # Define the path where to save the results.
    res_path = os.path.join(
        "/users/csb/huizing/Documents/PhD/Code/",
        "mowgli_reproducibility/evaluate/scores_seurat.pkl",
    )

    # Save the results.
    with open(res_path, "wb") as f:
        pickle.dump(scores_dict, f)

    console.log("Saved all of this! [bold green]")
