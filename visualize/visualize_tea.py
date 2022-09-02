# Imports
import matplotlib.pyplot as plt
import numpy as np
import pickle
import seaborn as sns
import os
import pandas as pd
import muon as mu
import scanpy as sc
import mofax
import anndata as ad

# Define the data folders.
data_folder = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/"
fig_folder = (
    "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/visualize/figures/"
)
w_folder = "/users/csb/huizing/Documents/PhD/Code/Mowgli/local_analysis/from_jz/w/"


# Load the data.
mdata = mu.read_h5mu(data_folder + "TEA/tea_preprocessed.h5mu.gz")

# Load the MOFA model.
mofa_model = mofax.mofa_model(data_folder + "TEA/tea_mofa_15.hdf5")
mdata.obsm["X_mofa"] = mofa_model.get_factors()

# Load the Mowgli model.
mdata.obsm["X_mowgli"] = np.load(
    w_folder + "tea_mowgli_cosine_50_0_05_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    allow_pickle=True,
).item()["W"]

# This is needed somehow.
mdata.uns = {}

# Make an object for the Mowgli embedding.
mowgli_embedding = ad.AnnData(mdata.obsm["X_mowgli"], obs=mdata.obs)
mowgli_embedding.obs[mdata["adt"].var_names] = mdata["adt"].X

# Make an object for the MOFA embedding.
mofa_embedding = ad.AnnData(mdata.obsm["X_mofa"], obs=mdata.obs)
mofa_embedding.obs[mdata["adt"].var_names] = mdata["adt"].X

# Compute neighbors for Mowgli and MOFA.
sc.pp.neighbors(mowgli_embedding, n_neighbors=25)
sc.pp.neighbors(mofa_embedding, n_neighbors=25)

# Compute UMAP for Mowgli and MOFA.
sc.tl.umap(mowgli_embedding)
sc.tl.umap(mofa_embedding)

# Compute Leiden for Mowgli and MOFA.
# The extra level of high resolution is useful to clean up the heatmap.
sc.tl.leiden(mowgli_embedding, resolution=0.2, key_added="leiden")
sc.tl.leiden(mowgli_embedding, resolution=3, key_added="leiden_precise")
sc.tl.leiden(mofa_embedding, resolution=0.2, key_added="leiden")
sc.tl.leiden(mofa_embedding, resolution=3, key_added="leiden_precise")

# Plot the unannotated UMAP plots with Leiden clusters.
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
sc.pl.umap(
    mowgli_embedding,
    color="leiden",
    alpha=0.7,
    legend_loc="on data",
    ax=axes[0],
    show=False,
)
sc.pl.umap(
    mofa_embedding,
    color="leiden",
    alpha=0.7,
    legend_loc="on data",
    ax=axes[1],
    show=False,
)
plt.savefig(fig_folder + "tea_leiden_unannotated.pdf")

# Annotate the MOFA emebdding.
cluster_names = {
    1: "B cells",
    8: "B cells",
    0: "CD4 T cells",
    2: "CD4 T cells",
    4: "CD8 T cells",
    6: "CD8 T cells",
    3: "Monocytes",
    5: "NK cells",
    7: "MAIT T cells",
    9: "Erythroid cells",
}
mofa_embedding.obs["leiden"] = [
    cluster_names[c] for c in mofa_embedding.obs["leiden"].cat.codes
]

# Annotate the Mowgli embedding.
cluster_names = {
    0: "B cells",
    1: "CD4 T cells",
    2: "CD4 T cells",
    3: "CD8 T cells",
    4: "CD4 T cells",
    5: "Monocytes",
    6: "NK cells",
    7: "MAIT T cells",
    8: "Erythroid cells",
}
mowgli_embedding.obs["leiden"] = [
    cluster_names[c] for c in mowgli_embedding.obs["leiden"].cat.codes
]

# Plot the annotated UMAP plots with Leiden clusters.
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
sc.pl.umap(
    mowgli_embedding,
    color="leiden",
    alpha=0.7,
    legend_loc="on data",
    ax=axes[0],
    show=False,
)
sc.pl.umap(
    mofa_embedding,
    color="leiden",
    alpha=0.7,
    legend_loc="on data",
    ax=axes[1],
    show=False,
)
plt.savefig(fig_folder + "tea_leiden_annotated.pdf")

# Compute dendrograms for leiden clusters (we use high resolution here!)
sc.tl.dendrogram(mowgli_embedding, groupby="leiden_precise")
sc.tl.dendrogram(mofa_embedding, groupby="leiden_precise")

# Rename the clusters based on their position in the dendrogram.
idx = mowgli_embedding.uns["dendrogram_leiden_precise"]["categories_idx_ordered"]
cats = mowgli_embedding.obs["leiden_precise"].cat.rename_categories(np.argsort(idx))
mowgli_embedding.obs["leiden_precise"] = cats

idx = mofa_embedding.uns["dendrogram_leiden_precise"]["categories_idx_ordered"]
cats = mofa_embedding.obs["leiden_precise"].cat.rename_categories(np.argsort(idx))
mofa_embedding.obs["leiden_precise"] = cats

# Rearrange the cells to mimic the precise dendrogram.
mowgli_embedding = mowgli_embedding[mowgli_embedding.obs["leiden_precise"].argsort()]
mofa_embedding = mofa_embedding[mofa_embedding.obs["leiden_precise"].argsort()]

# Compute the dendrogram for the less precise, annotated clusters.
sc.tl.dendrogram(mowgli_embedding, groupby="leiden")
sc.tl.dendrogram(mofa_embedding, groupby="leiden")

# Make heatmaps
sc.pl.heatmap(
    mowgli_embedding,
    var_names=mowgli_embedding.var_names[
        mowgli_embedding.X.std(0).argsort()[::-1][:15]
    ],
    groupby="leiden",
    dendrogram=True,
    show=False,
)
plt.savefig(fig_folder + "mowgli_tea_heatmap_annotated.pdf")

sc.pl.heatmap(
    mofa_embedding,
    var_names=mofa_embedding.var_names,
    groupby="leiden",
    dendrogram=True,
    show=False,
    vmin=-10,
    vmax=10,
    cmap="RdBu",
)
plt.savefig(fig_folder + "mofa_tea_heatmap_annotated.pdf")

# Plot UMAPS for all proteins.
sc.pl.umap(mowgli_embedding, color=mdata["adt"].var_names, show=False)
plt.savefig(fig_folder + "mowgli_tea_umap_proteins.pdf")

sc.pl.umap(mofa_embedding, color=mdata["adt"].var_names, show=False)
plt.savefig(fig_folder + "mofa_tea_umap_proteins.pdf")

# Plot UMAPS for all dimensions.
sc.pl.umap(
    mowgli_embedding,
    color=mowgli_embedding.var_names[mowgli_embedding.X.std(0).argsort()[::-1]],
    show=False,
    vmax=1,
)
plt.savefig(fig_folder + "mowgli_tea_umap_dims.pdf")

sc.pl.umap(
    mofa_embedding,
    color=mofa_embedding.var_names,
    show=False,
    vmin=-10,
    vmax=10,
    cmap="RdBu",
)
plt.savefig(fig_folder + "mofa_tea_umap_dims.pdf")
