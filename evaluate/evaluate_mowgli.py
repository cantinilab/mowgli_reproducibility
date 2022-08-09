######################################## IMPORTS #########################################

# Load libraries.
import scanpy as sc
import muon as mu
from rich.console import Console
import os
import numpy as np
from mowgli import score
import pickle

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


################################# EVALUATING MOWGLI ########################################

with console.status("[bold green]Evaluating Mowgli...") as status:

    # Intialize a dictionary for the scores.
    scores_dict = {}

    # Set the range of nearest neighbors.
    k_range = list(range(1, 30))

    # Set the range of resulotions.
    res_range = list(np.arange(0.1, 2, 0.1))

    for xp_name in ["pbmc_mowgli_50"]:

        # Initialise scores for this experiment.
        scores_dict[xp_name] = {}

        # Log the experiment name.
        console.log(f"Starting to compute scores for {xp_name} [bold green]")

        # Load the mowgli embedding.
        mowgli_path = os.path.join(
            "/users/csb/huizing/Documents/PhD/Code/",
            f"mowgli_reproducibility/data/10X_PBMC_10k/{xp_name}.npy",
        )
        mdata.obsm["X_mowgli"] = np.load(mowgli_path)
        mdata.uns = {}

        console.log("Loaded the embedding [bold green]")

        # Compute the silhouette score.
        scores_dict[xp_name]["Silhouette score"] = score.embedding_silhouette_score(
            mdata.obsm["X_mowgli"],
            mdata.obs["rna:celltype"],
        )

        console.log("Computed the silhouette score [bold green]")

        # Compute the kNN from the embedding.
        knn = score.embedding_to_knn(mdata.obsm["X_mowgli"], k=k_range[-1])

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
        sc.pp.neighbors(mdata, use_rep="X_mowgli")

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

    # Define the path where to save the results.
    res_path = os.path.join(
        "/users/csb/huizing/Documents/PhD/Code/",
        "mowgli_reproducibility/data/10X_PBMC_10k/scores_mowgli.pkl",
    )

    # Save the results.
    with open(res_path, "wb") as f:
        pickle.dump(scores_dict, f)

    console.log("Saved all of this! [bold green]")
