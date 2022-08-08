######################################## IMPORTS #########################################

# Load libraries.
import scanpy as sc
import muon as mu
from rich.console import Console
import os
import numpy as np
from mowgli import score
import mofax

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

################################# EVALUATING MOFA ########################################

with console.status("[bold green]Evaluating MOFA...") as status:

    # Intialize a dictionary for the scores.
    scores_dict = {}

    # Set the range of nearest neighbors.
    k_range = list(range(1, 30))

    # Set the range of resulotions.
    res_range = list(range(.1, 2, .1))

    for xp_name in ['pbmc_mofa_15', 'pbmc_mofa_30', 'pbmc_mofa_50']:

        # Initialise scores for this experiment.
        scores_dict[xp_name] = {}

        # Log the experiment name.
        console.log(f"Starting to compute scores for {xp_name} [bold green]")

        # Load the mofa embedding. TODO: mofax
        mofa_path = os.path.join(
            "/users/csb/huizing/Documents/PhD/Code/",
            f"mowgli_reproducibility/data/10X_PBMC_10k/{xp_name}.hdf5"
        )
        mofa_object = mofax.model_blabla(mofa_path)

        console.log("Loaded the model [bold green]")

        # TODO: put in obsm
        mdata.obsm["X_mofa"] = mofa_object.blablabla
        mdata.uns = {}

        # Compute the silhouette score.
        scores_dict[xp_name]['Silhouette score'] = score.embedding_silhouette_score(
            mdata.obsm["X_mofa"],
            mdata.obs["rna:celltype"],
        )

        console.log("Computed the silhouette score [bold green]")

        # Compute the purity score for varying k nearest neighbors.
        purity_scores = []
        for k in k_range:
            knn = score.embedding_to_knn(mdata.obsm["X_mofa"], k=k)
            purity_scores.append(score.knn_purity_score(knn, mdata.obs["rna:celltype"]))
        scores_dict[xp_name]['Purity scores'] = purity_scores
        scores_dict[xp_name]['k range'] = k_range

        console.log("Computed the purity scores. Phew! [bold green]")

        # Let Scanpy compute the kNN graph.
        sc.pp.neighbors(mdata, use_rep="X_mofa")

        # Compute the Leiden clustering and ARI for varying resolution.
        aris = []
        for res in res_range:
            sc.tl.leiden(mdata, res=res)
            aris.append(score.ARI(mdata.obs["rna:celltype"], mdata.obs["leiden"]))
        scores_dict[xp_name]['ARIs'] = aris
        scores_dict[xp_name]['res_range'] = res_range

        console.log("Computed the ARIs. Phew! [bold green]")

        # Try jaccard smoothing.
        jaccard_denoising(mdata)

        sc.tl.leiden(mdata)
        jaccard_aris = []
        for res in res_range:
            sc.tl.leiden(mdata, res=res)
            jaccard_aris.append(score.ARI(mdata.obs["rna:celltype"], mdata.obs["leiden"]))
        scores_dict[xp_name]['ARIs after denoising'] = jaccard_aris

        console.log("Computed the ARIs after denoising. Phew! [bold green]")

    # Define the path where to save the results.
    res_path = os.path.join(
        "/users/csb/huizing/Documents/PhD/Code/",
        "mowgli_reproducibility/data/10X_PBMC_10k/scores_mofa.pkl"
    )

    # Save the results.
    with open(res_path, 'wb') as f:
        pickle.dump(scores_dict, f)

    console.log("Saved all of this! [bold green]")