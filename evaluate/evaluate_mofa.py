######################################## IMPORTS #########################################

# Load libraries.
import scanpy as sc
import muon as mu
from rich.console import Console
import os
import numpy as np
from mowgli import score
import mofax
import pickle

console = Console()


# Define the data and figure folder.
data_folder = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/"

# Define the path where to save the results.
res_path = os.path.join(
    "/users/csb/huizing/Documents/PhD/Code/",
    "mowgli_reproducibility/evaluate/scores_mofa.pkl",
)

# Define data paths for different datasets.
data_path = {
    # "bmcite_mofa_15": data_folder + "BMCITE/bmcite_preprocessed.h5mu.gz",
    # "bmcite_mofa_30": data_folder + "BMCITE/bmcite_preprocessed.h5mu.gz",
    # "bmcite_mofa_50": data_folder + "BMCITE/bmcite_preprocessed.h5mu.gz",
    # "liu_mofa_5": data_folder + "Liu/liu_preprocessed.h5mu.gz",
    # "liu_mofa_15": data_folder + "Liu/liu_preprocessed.h5mu.gz",
    # "liu_mofa_30": data_folder + "Liu/liu_preprocessed.h5mu.gz",
    # "liu_mofa_50": data_folder + "Liu/liu_preprocessed.h5mu.gz",
    "sim1_mofa_5": data_folder + "Liu/liu_simulated_1.h5mu.gz",
    "sim1_mofa_15": data_folder + "Liu/liu_simulated_1.h5mu.gz",
    "sim1_mofa_30": data_folder + "Liu/liu_simulated_1.h5mu.gz",
    "sim1_mofa_50": data_folder + "Liu/liu_simulated_1.h5mu.gz",
    "sim2_mofa_5": data_folder + "Liu/liu_simulated_2.h5mu.gz",
    "sim2_mofa_15": data_folder + "Liu/liu_simulated_2.h5mu.gz",
    "sim2_mofa_30": data_folder + "Liu/liu_simulated_2.h5mu.gz",
    "sim2_mofa_50": data_folder + "Liu/liu_simulated_2.h5mu.gz",
    "sim3_mofa_5": data_folder + "Liu/liu_simulated_3.h5mu.gz",
    "sim3_mofa_15": data_folder + "Liu/liu_simulated_3.h5mu.gz",
    "sim3_mofa_30": data_folder + "Liu/liu_simulated_3.h5mu.gz",
    "sim3_mofa_50": data_folder + "Liu/liu_simulated_3.h5mu.gz",
    "sim4_mofa_5": data_folder + "Liu/liu_simulated_4.h5mu.gz",
    "sim4_mofa_15": data_folder + "Liu/liu_simulated_4.h5mu.gz",
    "sim4_mofa_30": data_folder + "Liu/liu_simulated_4.h5mu.gz",
    "sim4_mofa_50": data_folder + "Liu/liu_simulated_4.h5mu.gz",
    "sim5_mofa_5": data_folder + "Liu/liu_simulated_5.h5mu.gz",
    "sim5_mofa_15": data_folder + "Liu/liu_simulated_5.h5mu.gz",
    "sim5_mofa_30": data_folder + "Liu/liu_simulated_5.h5mu.gz",
    "sim5_mofa_50": data_folder + "Liu/liu_simulated_5.h5mu.gz",
    "sim6_mofa_5": data_folder + "Liu/liu_simulated_6.h5mu.gz",
    "sim6_mofa_15": data_folder + "Liu/liu_simulated_6.h5mu.gz",
    "sim6_mofa_30": data_folder + "Liu/liu_simulated_6.h5mu.gz",
    "sim6_mofa_50": data_folder + "Liu/liu_simulated_6.h5mu.gz",
    # "opcite_mofa_15": data_folder + "OPCITE/opcite_preprocessed.h5mu.gz",
    # "opcite_mofa_30": data_folder + "OPCITE/opcite_preprocessed.h5mu.gz",
    # "opcite_mofa_50": data_folder + "OPCITE/opcite_preprocessed.h5mu.gz",
    # "opmultiome_mofa_15": data_folder + "OP_multiome/opmultiome_preprocessed.h5mu.gz",
    # "opmultiome_mofa_30": data_folder + "OP_multiome/opmultiome_preprocessed.h5mu.gz",
    # "opmultiome_mofa_50": data_folder + "OP_multiome/opmultiome_preprocessed.h5mu.gz",
    # "pbmc_mofa_15": data_folder + "10X_PBMC_10k/pbmc_preprocessed.h5mu.gz",
    # "pbmc_mofa_30": data_folder + "10X_PBMC_10k/pbmc_preprocessed.h5mu.gz",
    # "pbmc_mofa_50": data_folder + "10X_PBMC_10k/pbmc_preprocessed.h5mu.gz",
    # "tea_mofa_15": data_folder + "TEA/tea_preprocessed.h5mu.gz",
    # "tea_mofa_30": data_folder + "TEA/tea_preprocessed.h5mu.gz",
    # "tea_mofa_50": data_folder + "TEA/tea_preprocessed.h5mu.gz",
}

mofa_path = {
    # "bmcite_mofa_15": data_folder + "BMCITE/bmcite_mofa_15.hdf5",
    # "bmcite_mofa_30": data_folder + "BMCITE/bmcite_mofa_30.hdf5",
    # "bmcite_mofa_50": data_folder + "BMCITE/bmcite_mofa_50.hdf5",
    # "liu_mofa_5": data_folder + "Liu/liu_mofa_5.hdf5",
    # "liu_mofa_15": data_folder + "Liu/liu_mofa_15.hdf5",
    # "liu_mofa_30": data_folder + "Liu/liu_mofa_30.hdf5",
    # "liu_mofa_50": data_folder + "Liu/liu_mofa_50.hdf5",
    "sim1_mofa_5": data_folder + "Liu/liu_simulated_1_mofa_5.hdf5",
    "sim1_mofa_15": data_folder + "Liu/liu_simulated_1_mofa_15.hdf5",
    "sim1_mofa_30": data_folder + "Liu/liu_simulated_1_mofa_30.hdf5",
    "sim1_mofa_50": data_folder + "Liu/liu_simulated_1_mofa_50.hdf5",
    "sim2_mofa_5": data_folder + "Liu/liu_simulated_2_mofa_5.hdf5",
    "sim2_mofa_15": data_folder + "Liu/liu_simulated_2_mofa_15.hdf5",
    "sim2_mofa_30": data_folder + "Liu/liu_simulated_2_mofa_30.hdf5",
    "sim2_mofa_50": data_folder + "Liu/liu_simulated_2_mofa_50.hdf5",
    "sim3_mofa_5": data_folder + "Liu/liu_simulated_3_mofa_5.hdf5",
    "sim3_mofa_15": data_folder + "Liu/liu_simulated_3_mofa_15.hdf5",
    "sim3_mofa_30": data_folder + "Liu/liu_simulated_3_mofa_30.hdf5",
    "sim3_mofa_50": data_folder + "Liu/liu_simulated_3_mofa_50.hdf5",
    "sim4_mofa_5": data_folder + "Liu/liu_simulated_4_mofa_5.hdf5",
    "sim4_mofa_15": data_folder + "Liu/liu_simulated_4_mofa_15.hdf5",
    "sim4_mofa_30": data_folder + "Liu/liu_simulated_4_mofa_30.hdf5",
    "sim4_mofa_50": data_folder + "Liu/liu_simulated_4_mofa_50.hdf5",
    "sim5_mofa_5": data_folder + "Liu/liu_simulated_5_mofa_5.hdf5",
    "sim5_mofa_15": data_folder + "Liu/liu_simulated_5_mofa_15.hdf5",
    "sim5_mofa_30": data_folder + "Liu/liu_simulated_5_mofa_30.hdf5",
    "sim5_mofa_50": data_folder + "Liu/liu_simulated_5_mofa_50.hdf5",
    "sim6_mofa_5": data_folder + "Liu/liu_simulated_6_mofa_5.hdf5",
    "sim6_mofa_15": data_folder + "Liu/liu_simulated_6_mofa_15.hdf5",
    "sim6_mofa_30": data_folder + "Liu/liu_simulated_6_mofa_30.hdf5",
    "sim6_mofa_50": data_folder + "Liu/liu_simulated_6_mofa_50.hdf5",
    # "opcite_mofa_15": data_folder + "OPCITE/opcite_mofa_15.hdf5",
    # "opcite_mofa_30": data_folder + "OPCITE/opcite_mofa_30.hdf5",
    # "opcite_mofa_50": data_folder + "OPCITE/opcite_mofa_50.hdf5",
    # "opmultiome_mofa_15": data_folder + "OP_multiome/opmultiome_mofa_15.hdf5",
    # "opmultiome_mofa_30": data_folder + "OP_multiome/opmultiome_mofa_30.hdf5",
    # "opmultiome_mofa_50": data_folder + "OP_multiome/opmultiome_mofa_50.hdf5",
    # "pbmc_mofa_15": data_folder + "10X_PBMC_10k/pbmc_mofa_15.hdf5",
    # "pbmc_mofa_30": data_folder + "10X_PBMC_10k/pbmc_mofa_30.hdf5",
    # "pbmc_mofa_50": data_folder + "10X_PBMC_10k/pbmc_mofa_50.hdf5",
    # "tea_mofa_15": data_folder + "TEA/tea_mofa_15.hdf5",
    # "tea_mofa_30": data_folder + "TEA/tea_mofa_30.hdf5",
    # "tea_mofa_50": data_folder + "TEA/tea_mofa_50.hdf5",
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


################################# EVALUATING MOFA ########################################

with console.status("[bold green]Evaluating MOFA...") as status:

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

        # Load the mofa embedding.
        mofa_object = mofax.mofa_model(mofa_path[xp_name])
        mdata.obsm["X_mofa"] = mofa_object.get_factors()
        mdata.uns = {}

        console.log("Loaded the model [bold green]")

        # Compute the silhouette score.
        scores_dict[xp_name]["Silhouette score"] = score.embedding_silhouette_score(
            mdata.obsm["X_mofa"],
            mdata.obs["rna:celltype"],
        )

        console.log("Computed the silhouette score [bold green]")

        # Compute the kNN from the embedding.
        knn = score.embedding_to_knn(mdata.obsm["X_mofa"], k=k_range[-1])

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
        sc.pp.neighbors(mdata, use_rep="X_mofa", n_neighbors=20)

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
