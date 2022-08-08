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

# TODO: different files for each method.

################################# EVALUATING MOFA ########################################

with console.status("[bold green]Evaluating MOFA...") as status:

    # TODO: init score dict.
    scores_dict = {}
    
    # TODO: faire varier les XP
    # Load the mofa embedding.
    mofa_path = os.path.join(
        "/users/csb/huizing/Documents/PhD/Code/",
        "mowgli_reproducibility/data/10X_PBMC_10k/pbmc_mofa_30.hdf5",
    )

    # TODO: mofax
    mofa_object = mofax.model_blabla

    # TODO: put in obsm
    mdata.obsm["X_mofa"] = np.load(nmf_path)
    mdata.uns = {}

    # Compute the silhouette score.
    sil = score.embedding_silhouette_score(
        mdata.obsm["X_mofa"],
        mdata.obs["rna:celltype"],
    )
    console.log("Silhouette score: [bold green]{:.2f}".format(sil))

    # Compute the purity score.
    knn = score.embedding_to_knn(mdata.obsm["X_mofa"])

    # TODO: faire varier k
    pur = score.knn_purity_score(knn, mdata.obs["rna:celltype"])
    console.log("Purity score: [bold green]{:.2f}".format(pur))

    # Compute the kNN graph.
    sc.pp.neighbors(mdata, use_rep="X_mofa")

    # Compute the Leiden clustering.
    # TODO: faire varier resolution
    sc.tl.leiden(mdata)

    # Compute the adjusted rand index.
    ari = score.ARI(mdata.obs["rna:celltype"], mdata.obs["leiden"])
    console.log("ARI: [bold green]{:.2f}".format(ari))

    # Try jaccard smoothing.
    jaccard_denoising(mdata)
    # TODO: faire varier resolution
    sc.tl.leiden(mdata)
    ari = score.ARI(mdata.obs["rna:celltype"], mdata.obs["leiden"])
    console.log("ARI after denoising: [bold green]{:.2f}".format(ari))

    # TODO: save somewhere. as mofa_results.pkl
    with open('saved_dictionary.pkl', 'wb') as f:
    pickle.dump(dictionary, f)        
    with open('saved_dictionary.pkl', 'rb') as f:
        loaded_dict = pickle.load(f)

################################## EVALUATING NMF ########################################

with console.status("[bold green]Evaluating NMF...") as status:

    # Load the NMF embedding.
    nmf_path = os.path.join(
        "/users/csb/huizing/Documents/PhD/Code/",
        "mowgli_reproducibility/data/10X_PBMC_10k/pbmc_nmf.npy",
    )
    mdata.obsm["X_nmf"] = np.load(nmf_path)
    mdata.uns = {}

    # Compute the silhouette score.
    sil = score.embedding_silhouette_score(
        mdata.obsm["X_nmf"],
        mdata.obs["rna:celltype"],
    )
    console.log("Silhouette score: [bold green]{:.2f}".format(sil))

    # Compute the purity score.
    knn = score.embedding_to_knn(mdata.obsm["X_nmf"])
    pur = score.knn_purity_score(knn, mdata.obs["rna:celltype"])
    console.log("Purity score: [bold green]{:.2f}".format(pur))

    # Compute the kNN graph.
    sc.pp.neighbors(mdata, use_rep="X_nmf")

    # Compute the Leiden clustering.
    sc.tl.leiden(mdata)

    # Compute the adjusted rand index.
    ari = score.ARI(mdata.obs["rna:celltype"], mdata.obs["leiden"])
    console.log("ARI: [bold green]{:.2f}".format(ari))

    # Try jaccard smoothing.
    jaccard_denoising(mdata)
    sc.tl.leiden(mdata)
    ari = score.ARI(mdata.obs["rna:celltype"], mdata.obs["leiden"])
    console.log("ARI after denoising: [bold green]{:.2f}".format(ari))
