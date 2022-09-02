# Imports
import matplotlib.pyplot as plt
import numpy as np
import pickle
import seaborn as sns
import os
import pandas as pd
from tueplots import axes as tue_axes
from tueplots import cycler as tue_cycler
from tueplots import fonts as tue_fonts
from tueplots.constants import markers as tue_markers
from tueplots.constants.color import palettes as tue_palettes

# Plots configuration.
plt.rcParams.update({"figure.dpi": 80})
plt.rcParams.update(tue_axes.spines(left=True, right=False, top=False, bottom=True))
plt.rcParams.update(tue_axes.grid())
plt.rcParams.update(tue_cycler.cycler(color=tue_palettes.high_contrast))
plt.rcParams.update(tue_axes.legend(shadow=False, frameon=False, fancybox=False))
plt.rcParams.update(tue_fonts.neurips2021_tex(family="sans-serif"))

# Define the paths where to read the results.
data_path = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/"
nmf_path = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/evaluate/scores_nmf.pkl"
mofa_path = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/evaluate/scores_mofa.pkl"
seurat_path = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/evaluate/scores_seurat.pkl"
mowgli_path = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/evaluate/scores_mowgli.pkl"

# Load the results from the pickle files.
with open(mofa_path, "rb") as f:
    scores_dict_mofa = pickle.load(f)
with open(nmf_path, "rb") as f:
    scores_dict_nmf = pickle.load(f)
with open(mowgli_path, "rb") as f:
    scores_dict_mowgli = pickle.load(f)
with open(seurat_path, "rb") as f:
    scores_dict_seurat = pickle.load(f)

# Fuse the dictionaries.
scores_dict = {
    **scores_dict_mofa,
    **scores_dict_nmf,
    **scores_dict_mowgli,
    **scores_dict_seurat,
}

# Turn the scores into a dataframe.
scores_df = pd.DataFrame(scores_dict).T

# Add a column for the method.
map_dict = {"mofa": "MOFA+", "nmf": "NMF", "mowgli": "Mowgli", "seurat": "Seurat"}
for x in map_dict:
    idx = scores_df.index.to_series().str.contains(x)
    scores_df.loc[idx, "Method"] = map_dict[x]

# Add a column for the dataset.
map_dict = {
    "pbmc": "10X PBMC",
    "liu": "Liu",
    "sim1": "Simulated 1",
    "sim2": "Simulated 2",
    "sim3": "Simulated 3",
    "sim4": "Simulated 4",
    "sim5": "Simulated 5",
    "sim6": "Simulated 6",
    "bmcite": "Bone Marrow CITE-seq",
    "opcite": "Open Problems CITE-seq",
    "opmultiome": "Open Problems Multiome",
}
for x in map_dict:
    idx = scores_df.index.to_series().str.contains(x)
    scores_df.loc[idx, "Dataset"] = map_dict[x]

# Add a column for the latent dimension.
# Careful, "5" is included in "50".
for x in ["5", "15", "30", "50"]:
    idx = scores_df.index.to_series().str.contains(x)
    scores_df.loc[idx, "Latent dimension"] = x

# Make a new dataframe with individual ARI scores for each resolution.
ari_res = []  # Initialize the list that will be turned into a dataframe.

# Iterate over experiments.
for xp_name in scores_df.index:

    # Iterate over resolutions.
    for i, res in enumerate(scores_df.loc[xp_name, "res_range"]):

        # Add the ARI to the list.
        ari_res.append(
            {
                "xp_name": xp_name,
                "Dataset": scores_df.loc[xp_name, "Dataset"],
                "Method": scores_df.loc[xp_name, "Method"],
                "Latent dimension": str(scores_df.loc[xp_name, "Latent dimension"]),
                "Resolution": res,
                "ARI": scores_df.loc[xp_name, "ARIs"][i],
            }
        )

# Turn the list into a dataframe.
ari_res = pd.DataFrame(ari_res)


# Make a new dataframe with individual purity scores depending on k.
purity_res = []  # Initialize the list that will be turned into a dataframe.

# Iterate over experiments.
for xp_name in scores_df.index:

    # Iterate over k nearest neighbours.
    for i, k in enumerate(scores_df.loc[xp_name, "k range"]):

        # Add the purity score to the list.
        purity_res.append(
            {
                "xp_name": xp_name,
                "Dataset": scores_df.loc[xp_name, "Dataset"],
                "Method": scores_df.loc[xp_name, "Method"],
                "Latent dimension": str(scores_df.loc[xp_name, "Latent dimension"]),
                "k": k,
                "Purity score": scores_df.loc[xp_name, "Purity scores"][i],
            }
        )

# Turn the list into a dataframe.
purity_res = pd.DataFrame(purity_res)


def plot_scores_dim_simulated(palette, method, path):

    # Define the datasets.
    datasets = [f"Simulated {i}" for i in range(1, 7)]

    # Define the order to plot datasets in.
    for i in range(len(datasets)):
        idx = scores_df["Dataset"].str.contains(datasets[i])
        scores_df.loc[idx, "Order"] = i

    # Define the subplots.
    fig = plt.figure(constrained_layout=True, figsize=(10, 10))
    axes = fig.subplot_mosaic(
        """
        ABC
        ADE
        AFG
        AHI
        AJK
        ALM
        """
    )

    # Visualize the silhouette score as a bar plot.
    idx = scores_df["Dataset"].str.contains("|".join(datasets))
    idx = idx & (scores_df["Method"] == method)
    sns.barplot(
        data=scores_df.loc[idx],
        y="Dataset",
        x="Silhouette score",
        hue="Latent dimension",
        hue_order=["5", "15", "30", "50"],
        palette=palette,
        order=datasets,
        ax=axes["A"],
    )
    axes["A"].legend(
        title="Latent dimension",
        bbox_to_anchor=(-0.1, 1),
        ncol=1,
    )
    axes["A"].set_yticklabels(axes["A"].get_yticklabels(), rotation=90, va="center")

    # Visualize the ARI as a line plot.
    ymin = ari_res["ARI"].min()
    ymax = ari_res["ARI"].max() + 0.05
    for i, ax in enumerate(
        [axes["B"], axes["D"], axes["F"], axes["H"], axes["J"], axes["L"]]
    ):
        idx = ari_res["Dataset"].str.contains(datasets[i])
        idx = idx & (ari_res["Method"] == method)
        sns.lineplot(
            data=ari_res.loc[idx],
            x="Resolution",
            y="ARI",
            hue="Latent dimension",
            hue_order=["5", "15", "30", "50"],
            palette=palette,
            ax=ax,
        )
        ax.get_legend().remove()
        ax.set(ylim=(ymin, ymax))
        if i < len(datasets) - 1:
            ax.set_xticklabels([])
            ax.set(xlabel=None)

    # Visualize the purity score as a line plot.
    ymin = purity_res["Purity score"].min()
    ymax = purity_res["Purity score"].max() + 0.05
    for i, ax in enumerate(
        [axes["C"], axes["E"], axes["G"], axes["I"], axes["K"], axes["M"]]
    ):
        idx = purity_res["Dataset"].str.contains(datasets[i])
        idx = idx & (purity_res["Method"] == method)
        sns.lineplot(
            data=purity_res.loc[idx],
            x="k",
            y="Purity score",
            hue="Latent dimension",
            hue_order=["5", "15", "30", "50"],
            palette=palette,
            ax=ax,
        )
        ax.get_legend().remove()
        ax.set(ylim=(0.45, ymax))
        if i < len(datasets) - 1:
            ax.set_xticklabels([])
            ax.set(xlabel=None)

    # Show a pretty grid.
    for i in axes:
        axes[i].grid()

    # Save the figure.
    fig.savefig(path)


def plot_scores_dim_real(palette, method, path):

    # Define the datasets.
    datasets = [
        "Liu",
        "10X PBMC",
        "Open Problems Multiome",
        "Open Problems CITE-seq",
        "Bone Marrow CITE-seq",
    ]

    # Define the order to plot datasets in.
    for i in range(len(datasets)):
        idx = scores_df["Dataset"].str.contains(datasets[i])
        scores_df.loc[idx, "Order"] = i

    # Define the subplots.
    fig = plt.figure(constrained_layout=True, figsize=(10, 10))
    axes = fig.subplot_mosaic(
        """
        ABC
        ADE
        AFG
        AHI
        AJK
        """
    )

    # Visualize the silhouette score as a bar plot.
    idx = scores_df["Dataset"].str.contains("|".join(datasets))
    idx = idx & (scores_df["Method"] == method)
    sns.barplot(
        data=scores_df.loc[idx],
        y="Dataset",
        x="Silhouette score",
        hue="Latent dimension",
        hue_order=["5", "15", "30", "50"],
        palette=palette,
        order=datasets,
        ax=axes["A"],
    )
    axes["A"].legend(
        title="Latent dimension",
        bbox_to_anchor=(-0.1, 1),
        ncol=1,
    )
    axes["A"].set_yticklabels(axes["A"].get_yticklabels(), rotation=90, va="center")

    # Visualize the ARI as a line plot.
    ymin = ari_res["ARI"].min()
    ymax = ari_res["ARI"].max() + 0.05
    for i, ax in enumerate(
        [axes["B"], axes["D"], axes["F"], axes["H"], axes["J"]]
    ):
        idx = ari_res["Dataset"].str.contains(datasets[i])
        idx = idx & (ari_res["Method"] == method)
        sns.lineplot(
            data=ari_res.loc[idx],
            x="Resolution",
            y="ARI",
            hue="Latent dimension",
            hue_order=["5", "15", "30", "50"],
            palette=palette,
            ax=ax,
        )
        ax.get_legend().remove()
        ax.set(ylim=(ymin, ymax))
        if i < len(datasets) - 1:
            ax.set_xticklabels([])
            ax.set(xlabel=None)

    # Visualize the purity score as a line plot.
    ymin = purity_res["Purity score"].min()
    ymax = purity_res["Purity score"].max() + 0.05
    for i, ax in enumerate(
        [axes["C"], axes["E"], axes["G"], axes["I"], axes["K"]]
    ):
        idx = purity_res["Dataset"].str.contains(datasets[i])
        idx = idx & (purity_res["Method"] == method)
        sns.lineplot(
            data=purity_res.loc[idx],
            x="k",
            y="Purity score",
            hue="Latent dimension",
            hue_order=["5", "15", "30", "50"],
            palette=palette,
            ax=ax,
        )
        ax.get_legend().remove()
        ax.set(ylim=(0.65, ymax))
        if i < len(datasets) - 1:
            ax.set_xticklabels([])
            ax.set(xlabel=None)

    # Show a pretty grid.
    for i in axes:
        axes[i].grid()

    # Save the figure.
    fig.savefig(path)


# Plot the scores for Mowgli for simulated data, per dimension.
path = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/visualize/figures/scores_dim_simulated_mowgli.pdf"
palette = sns.color_palette("flare", n_colors=4)
plot_scores_dim_simulated(palette, "Mowgli", path)

# Plot the scores for MOFA+ for simulated data, per dimension.
path = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/visualize/figures/scores_dim_simulated_mofa.pdf"
palette = sns.color_palette("flare", n_colors=4)
plot_scores_dim_simulated(palette, "MOFA+", path)

# Plot the scores for NMF for simulated data, per dimension.
path = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/visualize/figures/scores_dim_simulated_nmf.pdf"
palette = sns.color_palette("flare", n_colors=4)
plot_scores_dim_simulated(palette, "NMF", path)

# Plot the scores for Mowgli for real data, per dimension.
path = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/visualize/figures/scores_dim_real_mowgli.pdf"
palette = sns.color_palette("flare", n_colors=4)
plot_scores_dim_real(palette, "Mowgli", path)

# Plot the scores for MOFA+ for real data, per dimension.
path = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/visualize/figures/scores_dim_real_mofa.pdf"
palette = sns.color_palette("flare", n_colors=4)
plot_scores_dim_real(palette, "MOFA+", path)

# Plot the scores for NMF for real data, per dimension.
path = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/visualize/figures/scores_dim_real_nmf.pdf"
palette = sns.color_palette("flare", n_colors=4)
plot_scores_dim_real(palette, "NMF", path)
