######################################## IMPORTS #########################################

# Load libraries.
import scanpy as sc
import muon as mu
from rich.console import Console
import os
import numpy as np
from mowgli import score
import pickle
import sys

# Define the data and figure folder.
data_folder = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/"

# Define the path where to save the results.
res_path = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/evaluate/scores_nmf.pkl"

# Define data paths for different datasets.
data_path = {
    # "liu_nmf_5": data_folder + "Liu/liu_preprocessed.h5mu.gz",
    # "sim1_nmf_5": data_folder + "Liu/liu_simulated_1.h5mu.gz",
    # "sim2_nmf_5": data_folder + "Liu/liu_simulated_2.h5mu.gz",
    # "sim3_nmf_5": data_folder + "Liu/liu_simulated_3.h5mu.gz",
    "sim4_nmf_5": data_folder + "Liu/liu_simulated_4.h5mu.gz",
    "sim5_nmf_5": data_folder + "Liu/liu_simulated_5.h5mu.gz",
    "sim6_nmf_5": data_folder + "Liu/liu_simulated_6.h5mu.gz",
    
    # "bmcite_nmf_15": data_folder + "BMCITE/bmcite_preprocessed.h5mu.gz",
    # "liu_nmf_15": data_folder + "Liu/liu_preprocessed.h5mu.gz",
    # "sim1_nmf_15": data_folder + "Liu/liu_simulated_1.h5mu.gz",
    # "sim2_nmf_15": data_folder + "Liu/liu_simulated_2.h5mu.gz",
    # "sim3_nmf_15": data_folder + "Liu/liu_simulated_3.h5mu.gz",
    "sim4_nmf_15": data_folder + "Liu/liu_simulated_4.h5mu.gz",
    "sim5_nmf_15": data_folder + "Liu/liu_simulated_5.h5mu.gz",
    "sim6_nmf_15": data_folder + "Liu/liu_simulated_6.h5mu.gz",
    # "opcite_nmf_15": data_folder + "OPCITE/opcite_preprocessed.h5mu.gz",
    # "opmultiome_nmf_15": data_folder + "OP_multiome/opmultiome_preprocessed.h5mu.gz",
    # "pbmc_nmf_15": data_folder + "10X_PBMC_10k/pbmc_preprocessed.h5mu.gz",
    # "tea_nmf_15": data_folder + "TEA/tea_preprocessed.h5mu.gz",

    # "bmcite_nmf_30": data_folder + "BMCITE/bmcite_preprocessed.h5mu.gz",
    # "liu_nmf_30": data_folder + "Liu/liu_preprocessed.h5mu.gz",
    # "sim1_nmf_30": data_folder + "Liu/liu_simulated_1.h5mu.gz",
    # "sim2_nmf_30": data_folder + "Liu/liu_simulated_2.h5mu.gz",
    # "sim3_nmf_30": data_folder + "Liu/liu_simulated_3.h5mu.gz",
    "sim4_nmf_30": data_folder + "Liu/liu_simulated_4.h5mu.gz",
    "sim5_nmf_30": data_folder + "Liu/liu_simulated_5.h5mu.gz",
    "sim6_nmf_30": data_folder + "Liu/liu_simulated_6.h5mu.gz",
    # "opcite_nmf_30": data_folder + "OPCITE/opcite_preprocessed.h5mu.gz",
    # "opmultiome_nmf_30": data_folder + "OP_multiome/opmultiome_preprocessed.h5mu.gz",
    # "pbmc_nmf_30": data_folder + "10X_PBMC_10k/pbmc_preprocessed.h5mu.gz",
    # "tea_nmf_30": data_folder + "TEA/tea_preprocessed.h5mu.gz",

    # "bmcite_nmf_50": data_folder + "BMCITE/bmcite_preprocessed.h5mu.gz",
    # "liu_nmf_50": data_folder + "Liu/liu_preprocessed.h5mu.gz",
    # "sim1_nmf_50": data_folder + "Liu/liu_simulated_1.h5mu.gz",
    # "sim2_nmf_50": data_folder + "Liu/liu_simulated_2.h5mu.gz",
    # "sim3_nmf_50": data_folder + "Liu/liu_simulated_3.h5mu.gz",
    "sim4_nmf_50": data_folder + "Liu/liu_simulated_4.h5mu.gz",
    "sim5_nmf_50": data_folder + "Liu/liu_simulated_5.h5mu.gz",
    "sim6_nmf_50": data_folder + "Liu/liu_simulated_6.h5mu.gz",
    # "opcite_nmf_50": data_folder + "OPCITE/opcite_preprocessed.h5mu.gz",
    # "opmultiome_nmf_50": data_folder + "OP_multiome/opmultiome_preprocessed.h5mu.gz",
    # "pbmc_nmf_50": data_folder + "10X_PBMC_10k/pbmc_preprocessed.h5mu.gz",
    # "tea_nmf_50": data_folder + "TEA/tea_preprocessed.h5mu.gz",
}

nmf_path = {
    # "liu_nmf_5": data_folder + "Liu/liu_nmf_5.npy",
    # "sim1_nmf_5": data_folder + "Liu/liu_simulated_1_nmf_5.npy",
    # "sim2_nmf_5": data_folder + "Liu/liu_simulated_2_nmf_5.npy",
    # "sim3_nmf_5": data_folder + "Liu/liu_simulated_3_nmf_5.npy",
    "sim4_nmf_5": data_folder + "Liu/liu_simulated_4_nmf_5.npy",
    "sim5_nmf_5": data_folder + "Liu/liu_simulated_5_nmf_5.npy",
    "sim6_nmf_5": data_folder + "Liu/liu_simulated_6_nmf_5.npy",
    
    # "bmcite_nmf_15": data_folder + "BMCITE/bmcite_nmf_15.npy",
    # "liu_nmf_15": data_folder + "Liu/liu_nmf_15.npy",
    # "sim1_nmf_15": data_folder + "Liu/liu_simulated_1_nmf_15.npy",
    # "sim2_nmf_15": data_folder + "Liu/liu_simulated_2_nmf_15.npy",
    # "sim3_nmf_15": data_folder + "Liu/liu_simulated_3_nmf_15.npy",
    "sim4_nmf_15": data_folder + "Liu/liu_simulated_4_nmf_15.npy",
    "sim5_nmf_15": data_folder + "Liu/liu_simulated_5_nmf_15.npy",
    "sim6_nmf_15": data_folder + "Liu/liu_simulated_6_nmf_15.npy",
    # "opcite_nmf_15": data_folder + "OPCITE/opcite_nmf_15.npy",
    # "opmultiome_nmf_15": data_folder + "OP_multiome/opmultiome_nmf_15.npy",
    # "pbmc_nmf_15": data_folder + "10X_PBMC_10k/pbmc_nmf_15.npy",
    # "tea_nmf_15": data_folder + "TEA/tea_nmf_15.npy",

    # "bmcite_nmf_30": data_folder + "BMCITE/bmcite_nmf_30.npy",
    # "liu_nmf_30": data_folder + "Liu/liu_nmf_30.npy",
    # "sim1_nmf_30": data_folder + "Liu/liu_simulated_1_nmf_30.npy",
    # "sim2_nmf_30": data_folder + "Liu/liu_simulated_2_nmf_30.npy",
    # "sim3_nmf_30": data_folder + "Liu/liu_simulated_3_nmf_30.npy",
    "sim4_nmf_30": data_folder + "Liu/liu_simulated_4_nmf_30.npy",
    "sim5_nmf_30": data_folder + "Liu/liu_simulated_5_nmf_30.npy",
    "sim6_nmf_30": data_folder + "Liu/liu_simulated_6_nmf_30.npy",
    # "opcite_nmf_30": data_folder + "OPCITE/opcite_nmf_30.npy",
    # "opmultiome_nmf_30": data_folder + "OP_multiome/opmultiome_nmf_30.npy",
    # "pbmc_nmf_30": data_folder + "10X_PBMC_10k/pbmc_nmf_30.npy",
    # "tea_nmf_30": data_folder + "TEA/tea_nmf_30.npy",

    # "bmcite_nmf_50": data_folder + "BMCITE/bmcite_nmf_50.npy",
    # "liu_nmf_50": data_folder + "Liu/liu_nmf_50.npy",
    # "sim1_nmf_50": data_folder + "Liu/liu_simulated_1_nmf_50.npy",
    # "sim2_nmf_50": data_folder + "Liu/liu_simulated_2_nmf_50.npy",
    # "sim3_nmf_50": data_folder + "Liu/liu_simulated_3_nmf_50.npy",
    "sim4_nmf_50": data_folder + "Liu/liu_simulated_4_nmf_50.npy",
    "sim5_nmf_50": data_folder + "Liu/liu_simulated_5_nmf_50.npy",
    "sim6_nmf_50": data_folder + "Liu/liu_simulated_6_nmf_50.npy",
    # "opcite_nmf_50": data_folder + "OPCITE/opcite_nmf_50.npy",
    # "opmultiome_nmf_50": data_folder + "OP_multiome/opmultiome_nmf_50.npy",
    # "pbmc_nmf_50": data_folder + "10X_PBMC_10k/pbmc_nmf_50.npy",
    # "tea_nmf_50": data_folder + "TEA/tea_nmf_50.npy",
}

console = Console()

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


################################# EVALUATING NMF ########################################

with console.status("[bold green]Evaluating NMF...") as status:

    # Intialize a dictionary for the scores.
    scores_dict = {}
    with open(res_path, "rb") as f:
        scores_dict = pickle.load(f)

    # Set the range of nearest neighbors.
    k_range = list(range(1, 30))

    # Set the range of resulotions.
    res_range = list(np.arange(0.1, 2, 0.1))

    previous_path = ""

    for xp_name in data_path:

        # Load the original data.
        console.log(f"Loading data for {xp_name} [bold green]")
        if previous_path != data_path[xp_name]:
            mdata = mu.read_h5mu(data_path[xp_name])
            previous_path = data_path[xp_name]
        console.log("Data loaded.")

        # Initialise scores for this experiment.
        scores_dict[xp_name] = {}

        # Log the experiment name.
        console.log(f"Starting to compute scores for {xp_name} [bold green]")

        # Load the nmf embedding.
        mdata.obsm["X_nmf"] = np.load(nmf_path[xp_name])
        mdata.uns = {}

        console.log("Loaded the embedding [bold green]")

        # Compute the silhouette score.
        scores_dict[xp_name]["Silhouette score"] = score.embedding_silhouette_score(
            mdata.obsm["X_nmf"],
            mdata.obs["rna:celltype"],
        )

        console.log("Computed the silhouette score [bold green]")

        # Compute the kNN from the embedding.
        knn = score.embedding_to_knn(mdata.obsm["X_nmf"], k=k_range[-1])

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

        # Let Scanpy compute the kNN graph.
        sc.pp.neighbors(mdata, use_rep="X_nmf", n_neighbors=20)

        # Compute the Leiden clustering and ARI for varying resolution.
        aris = []
        for res in res_range:

            # Log the value of resolution.
            console.log(f"Computing ARI for resolution={res} [bold green]")

            # Compute the ARI.
            sc.tl.leiden(mdata, resolution=res)
            aris.append(score.ARI(mdata.obs["rna:celltype"], mdata.obs["leiden"]))

        scores_dict[xp_name]["ARIs"] = aris
        scores_dict[xp_name]["res_range"] = res_range

        console.log("Computed the ARIs. Phew! [bold green]")

        # Try jaccard smoothing.
        jaccard_denoising(mdata)

        sc.tl.leiden(mdata)
        jaccard_aris = []
        for res in res_range:

            # Log the value of resolution.
            console.log(f"Computing ARI for resolution={res} [bold green]")

            # Compute the ARI.
            sc.tl.leiden(mdata, resolution=res)
            s = score.ARI(mdata.obs["rna:celltype"], mdata.obs["leiden"])
            jaccard_aris.append(s)

        scores_dict[xp_name]["ARIs after denoising"] = jaccard_aris

        console.log("Computed the ARIs after denoising. Phew! [bold green]")

    # Save the results.
    with open(res_path, "wb") as f:
        pickle.dump(scores_dict, f)

    console.log("Saved all of this! [bold green]")
