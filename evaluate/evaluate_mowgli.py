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
res_folder = "/users/csb/huizing/Documents/PhD/Code/Mowgli/local_analysis/from_jz/w/"

# Define data paths for different datasets.
data_path = {
    "liu_mowgli_cosine_5": data_folder + "Liu/liu_preprocessed.h5mu.gz",
    "liu_mowgli_cosine_15": data_folder + "Liu/liu_preprocessed.h5mu.gz",
    "liu_mowgli_cosine_30": data_folder + "Liu/liu_preprocessed.h5mu.gz",
    "liu_mowgli_cosine_50": data_folder + "Liu/liu_preprocessed.h5mu.gz",
    
    "sim1_mowgli_cosine_5": data_folder + "Liu/liu_simulated_1.h5mu.gz",
    "sim1_mowgli_cosine_15": data_folder + "Liu/liu_simulated_1.h5mu.gz",
    "sim1_mowgli_cosine_30": data_folder + "Liu/liu_simulated_1.h5mu.gz",
    "sim1_mowgli_cosine_50": data_folder + "Liu/liu_simulated_1.h5mu.gz",

    "sim2_mowgli_cosine_5": data_folder + "Liu/liu_simulated_2.h5mu.gz",
    "sim2_mowgli_cosine_15": data_folder + "Liu/liu_simulated_2.h5mu.gz",
    "sim2_mowgli_cosine_30": data_folder + "Liu/liu_simulated_2.h5mu.gz",
    "sim2_mowgli_cosine_50": data_folder + "Liu/liu_simulated_2.h5mu.gz",

    "sim3_mowgli_cosine_5": data_folder + "Liu/liu_simulated_3.h5mu.gz",
    "sim3_mowgli_cosine_15": data_folder + "Liu/liu_simulated_3.h5mu.gz",
    "sim3_mowgli_cosine_30": data_folder + "Liu/liu_simulated_3.h5mu.gz",
    "sim3_mowgli_cosine_50": data_folder + "Liu/liu_simulated_3.h5mu.gz",

    "sim4_mowgli_cosine_5": data_folder + "Liu/liu_simulated_4.h5mu.gz",
    "sim4_mowgli_cosine_15": data_folder + "Liu/liu_simulated_4.h5mu.gz",
    "sim4_mowgli_cosine_30": data_folder + "Liu/liu_simulated_4.h5mu.gz",
    "sim4_mowgli_cosine_50": data_folder + "Liu/liu_simulated_4.h5mu.gz",

    "sim5_mowgli_cosine_5": data_folder + "Liu/liu_simulated_5.h5mu.gz",
    "sim5_mowgli_cosine_15": data_folder + "Liu/liu_simulated_5.h5mu.gz",
    "sim5_mowgli_cosine_30": data_folder + "Liu/liu_simulated_5.h5mu.gz",
    "sim5_mowgli_cosine_50": data_folder + "Liu/liu_simulated_5.h5mu.gz",

    "sim6_mowgli_cosine_5": data_folder + "Liu/liu_simulated_6.h5mu.gz",
    "sim6_mowgli_cosine_15": data_folder + "Liu/liu_simulated_6.h5mu.gz",
    "sim6_mowgli_cosine_30": data_folder + "Liu/liu_simulated_6.h5mu.gz",
    "sim6_mowgli_cosine_50": data_folder + "Liu/liu_simulated_6.h5mu.gz",

    # "pbmc_mowgli_cosine_15": data_folder + "10X_PBMC_10k/pbmc_preprocessed.h5mu.gz",
    # "pbmc_mowgli_cosine_30": data_folder + "10X_PBMC_10k/pbmc_preprocessed.h5mu.gz",
    # "pbmc_mowgli_cosine_50": data_folder + "10X_PBMC_10k/pbmc_preprocessed.h5mu.gz",

    # "opmultiome_mowgli_cosine_15": data_folder + "OP_multiome/opmultiome_preprocessed.h5mu.gz",
    # "opmultiome_mowgli_cosine_30": data_folder + "OP_multiome/opmultiome_preprocessed.h5mu.gz",
    # "opmultiome_mowgli_cosine_50": data_folder + "OP_multiome/opmultiome_preprocessed.h5mu.gz",

    # "opcite_mowgli_cosine_15": data_folder + "OPCITE/opcite_preprocessed.h5mu.gz",
    # "opcite_mowgli_cosine_30": data_folder + "OPCITE/opcite_preprocessed.h5mu.gz",
    # "opcite_mowgli_cosine_50": data_folder + "OPCITE/opcite_preprocessed.h5mu.gz",

    # "bmcite_mowgli_cosine_15": data_folder + "BMCITE/bmcite_preprocessed.h5mu.gz",
    # "bmcite_mowgli_cosine_30": data_folder + "BMCITE/bmcite_preprocessed.h5mu.gz",
    # "bmcite_mowgli_cosine_50": data_folder + "BMCITE/bmcite_preprocessed.h5mu.gz",
}

mowgli_path = {
    "liu_mowgli_cosine_5": res_folder + "liu_mowgli_cosine_5_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    "liu_mowgli_cosine_15": res_folder + "liu_mowgli_cosine_15_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    "liu_mowgli_cosine_30": res_folder + "liu_mowgli_cosine_30_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    "liu_mowgli_cosine_50": res_folder + "liu_mowgli_cosine_50_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    
    "sim1_mowgli_cosine_5": res_folder + "liu_sim_1_mowgli_cosine_5_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    "sim1_mowgli_cosine_15": res_folder + "liu_sim_1_mowgli_cosine_15_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    "sim1_mowgli_cosine_30": res_folder + "liu_sim_1_mowgli_cosine_30_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    "sim1_mowgli_cosine_50": res_folder + "liu_sim_1_mowgli_cosine_50_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",

    "sim2_mowgli_cosine_5": res_folder + "liu_sim_2_mowgli_cosine_5_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    "sim2_mowgli_cosine_15": res_folder + "liu_sim_2_mowgli_cosine_15_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    "sim2_mowgli_cosine_30": res_folder + "liu_sim_2_mowgli_cosine_30_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    "sim2_mowgli_cosine_50": res_folder + "liu_sim_2_mowgli_cosine_50_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",

    "sim3_mowgli_cosine_5": res_folder + "liu_sim_3_mowgli_cosine_5_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    "sim3_mowgli_cosine_15": res_folder + "liu_sim_3_mowgli_cosine_15_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    "sim3_mowgli_cosine_30": res_folder + "liu_sim_3_mowgli_cosine_30_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    "sim3_mowgli_cosine_50": res_folder + "liu_sim_3_mowgli_cosine_50_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",

    "sim4_mowgli_cosine_5": res_folder + "liu_sim_4_mowgli_cosine_5_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    "sim4_mowgli_cosine_15": res_folder + "liu_sim_4_mowgli_cosine_15_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    "sim4_mowgli_cosine_30": res_folder + "liu_sim_4_mowgli_cosine_30_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    "sim4_mowgli_cosine_50": res_folder + "liu_sim_4_mowgli_cosine_50_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",

    "sim5_mowgli_cosine_5": res_folder + "liu_sim_5_mowgli_cosine_5_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    "sim5_mowgli_cosine_15": res_folder + "liu_sim_5_mowgli_cosine_15_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    "sim5_mowgli_cosine_30": res_folder + "liu_sim_5_mowgli_cosine_30_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    "sim5_mowgli_cosine_50": res_folder + "liu_sim_5_mowgli_cosine_50_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",

    "sim6_mowgli_cosine_5": res_folder + "liu_sim_6_mowgli_cosine_5_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    "sim6_mowgli_cosine_15": res_folder + "liu_sim_6_mowgli_cosine_15_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    "sim6_mowgli_cosine_30": res_folder + "liu_sim_6_mowgli_cosine_30_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    "sim6_mowgli_cosine_50": res_folder + "liu_sim_6_mowgli_cosine_50_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy",

    # "pbmc_mowgli_cosine_15": res_folder + "pbmc_mowgli_cosine_15_0_05_rna_0_01_atac_0_1_adt_0_001_0_001.npy",
    # "pbmc_mowgli_cosine_30": res_folder + "pbmc_mowgli_cosine_30_0_05_rna_0_01_atac_0_1_adt_0_001_0_001.npy",
    # "pbmc_mowgli_cosine_50": res_folder + "pbmc_mowgli_cosine_50_0_05_rna_0_01_atac_0_1_adt_0_001_0_001.npy",

    # "opmultiome_mowgli_cosine_15": res_folder + "opmultiome_mowgli_cosine_15_0_05_rna_0_01_atac_0_1_adt_0_001_0_001.npy",
    # "opmultiome_mowgli_cosine_30": res_folder + "opmultiome_mowgli_cosine_30_0_05_rna_0_01_atac_0_1_adt_0_001_0_001.npy",
    # "opmultiome_mowgli_cosine_50": res_folder + "opmultiome_mowgli_cosine_50_0_05_rna_0_01_atac_0_1_adt_0_001_0_001.npy",

    # "opcite_mowgli_cosine_15": res_folder + "opcite_mowgli_cosine_15_0_05_0_01_0_001.npy",
    # "opcite_mowgli_cosine_30": res_folder + "opcite_mowgli_cosine_30_0_05_0_01_0_001.npy",
    # "opcite_mowgli_cosine_50": res_folder + "opcite_mowgli_cosine_50_0_05_0_01_0_001.npy",

    # "bmcite_mowgli_cosine_15": res_folder + "bmcite_mowgli_cosine_15_0_05_0_01_0_001.npy",
    # "bmcite_mowgli_cosine_30": res_folder + "bmcite_mowgli_cosine_30_0_05_0_01_0_001.npy",
    # "bmcite_mowgli_cosine_50": res_folder + "bmcite_mowgli_cosine_50_0_05_0_01_0_001.npy",
}

# Define the path where to save the results.
res_path = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/evaluate/scores_mowgli.pkl"

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

            if d < 1/15: # Pruning threshold
                d = 0

            # Reweight the kNN.
            adata.obsp["connectivities"][i, j] = d
            adata.obsp["connectivities"][j, i] = d


################################# EVALUATING mowgli ########################################

with console.status("[bold green]Evaluating mowgli...") as status:

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

        # Load the mowgli embedding.
        mdata.obsm["X_mowgli"] = np.load(
            mowgli_path[xp_name], allow_pickle=True
        ).item()["W"]
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
        sc.pp.neighbors(mdata, use_rep="X_mowgli", n_neighbors=20)

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
