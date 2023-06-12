import hydra
from omegaconf import DictConfig


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def my_app(cfg: DictConfig) -> None:

    # Imports
    import pickle

    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns
    from tueplots import axes as tue_axes
    from tueplots import cycler as tue_cycler
    from tueplots import fonts as tue_fonts
    from tueplots.constants.color import palettes as tue_palettes

    # Plots configuration.
    plt.rcParams.update({"figure.dpi": 80})
    plt.rcParams.update(tue_axes.spines(left=True, right=False, top=False, bottom=True))
    plt.rcParams.update(tue_axes.grid())
    plt.rcParams.update(tue_cycler.cycler(color=tue_palettes.high_contrast))
    plt.rcParams.update(tue_axes.legend(shadow=False, frameon=False, fancybox=False))
    plt.rcParams.update(tue_fonts.neurips2021_tex(family="sans-serif"))

    # Define the paths where to read the results.
    nmf_path = f"{cfg.evaluate_path}scores_nmf.pkl"
    mofa_path = f"{cfg.evaluate_path}scores_mofa.pkl"
    seurat_path = f"{cfg.evaluate_path}scores_seurat.pkl"
    mowgli_path = f"{cfg.evaluate_path}scores_mowgli.pkl"
    cobolt_path = f"{cfg.evaluate_path}scores_cobolt.pkl"
    multigrate_path = f"{cfg.evaluate_path}scores_multigrate.pkl"

    # Load the results from the pickle files.
    with open(mofa_path, "rb") as f:
        scores_dict_mofa = pickle.load(f)
    with open(nmf_path, "rb") as f:
        scores_dict_nmf = pickle.load(f)
    with open(mowgli_path, "rb") as f:
        scores_dict_mowgli = pickle.load(f)
    with open(seurat_path, "rb") as f:
        scores_dict_seurat = pickle.load(f)
    with open(cobolt_path, "rb") as f:
        scores_dict_cobolt = pickle.load(f)
    with open(multigrate_path, "rb") as f:
        scores_dict_multigrate = pickle.load(f)

    # Fuse the dictionaries.
    scores_dict = {
        **scores_dict_mofa,
        **scores_dict_nmf,
        **scores_dict_mowgli,
        **scores_dict_seurat,
        **scores_dict_cobolt,
        **scores_dict_multigrate,
    }

    # Turn the scores into a dataframe.
    scores_df = pd.DataFrame(scores_dict).T

    # Add a column for the method.
    map_dict = {
        "mofa": "MOFA+",
        "nmf": "NMF",
        "mowgli": "Mowgli",
        "seurat": "Seurat",
        "cobolt": "Cobolt",
        "multigrate": "Multigrate",
    }
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
        "liu_simulated_1": "Simulated 1",
        "liu_simulated_2": "Simulated 2",
        "liu_simulated_3": "Simulated 3",
        "liu_simulated_4": "Simulated 4",
        "liu_simulated_5": "Simulated 5",
        "liu_simulated_6": "Simulated 6",
        "bmcite": "Bone Marrow CITE-seq",
        "opcite": "Open Problems CITE-seq",
        "opmultiome": "Open Problems Multiome",
    }
    for x in map_dict:
        idx = scores_df.index.to_series().str.contains(x)
        scores_df.loc[idx, "Dataset"] = map_dict[x]

    # Add a column for the latent dimension.
    # Careful, "5" is included in "50", but since 50 come later, this is corrected.
    for x in ["5", "15", "30", "50", "100"]:
        idx = scores_df.index.to_series().str.contains(x)
        scores_df.loc[idx, "Latent dimension"] = x

    # Add a column for the H_rna regularization.
    for x in [1.0, 0.1, 0.01, 0.001]:
        idx = scores_df.index.to_series().str.contains(f"rna_{x}".replace(".", "_"))
        scores_df.loc[idx, "H_rna regularization"] = x

    # Add a column for the H_atac regularization.
    for x in [1.0, 0.1, 0.01, 0.001]:
        idx = scores_df.index.to_series().str.contains(f"atac_{x}".replace(".", "_"))
        scores_df.loc[idx, "H_atac regularization"] = x

    # Add a column for the H_adt regularization.
    for x in [1.0, 0.1, 0.01, 0.001]:
        idx = scores_df.index.to_series().str.contains(f"adt_{x}".replace(".", "_"))
        scores_df.loc[idx, "H_adt regularization"] = x

    # Add a column for the w_regularization.
    for x in [0.01, 0.001, 0.0001]:
        idx = scores_df.index.to_series().str.endswith(f"_{x}".replace(".", "_"))
        scores_df.loc[idx, "w_regularization"] = x

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
                    "H_rna regularization": scores_df.loc[
                        xp_name, "H_rna regularization"
                    ],
                    "H_atac regularization": scores_df.loc[
                        xp_name, "H_atac regularization"
                    ],
                    "H_adt regularization": scores_df.loc[
                        xp_name, "H_adt regularization"
                    ],
                    "w_regularization": scores_df.loc[xp_name, "w_regularization"],
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
                    "H_rna regularization": scores_df.loc[
                        xp_name, "H_rna regularization"
                    ],
                    "H_atac regularization": scores_df.loc[
                        xp_name, "H_atac regularization"
                    ],
                    "H_adt regularization": scores_df.loc[
                        xp_name, "H_adt regularization"
                    ],
                    "w_regularization": scores_df.loc[xp_name, "w_regularization"],
                }
            )

    # Turn the list into a dataframe.
    purity_res = pd.DataFrame(purity_res)

    # For each dataset and method, get the latent dimension of the best silhouette score.
    best_latent_dim_sil = []
    for dataset in scores_df["Dataset"].unique():
        for method in scores_df["Method"].unique():
            idx = (scores_df["Dataset"] == dataset) & (scores_df["Method"] == method)
            if method == "Mowgli":
                idx &= scores_df["H_rna regularization"] == 0.01
                idx &= scores_df["H_atac regularization"] == 0.1
                idx &= scores_df["H_adt regularization"] == 0.01
                idx &= scores_df["w_regularization"] == 0.001
            best_latent_dim_sil.append(
                {
                    "Dataset": dataset,
                    "Method": method,
                    "Latent dimension": scores_df.loc[idx, "Latent dimension"][
                        scores_df.loc[idx, "Silhouette score"].astype("float").argmax()
                    ],
                    "Silhouette score": scores_df.loc[idx, "Silhouette score"].max(),
                }
            )
    best_latent_dim_sil = pd.DataFrame(best_latent_dim_sil)

    # Get the maximum ARI for each experiment.
    xp_max_ari = ari_res.groupby("xp_name")["ARI"].max()

    # Add this as a column to the dataframe.
    ari_res["max_ARI"] = ari_res["xp_name"].map(xp_max_ari)

    # For each dataset and method, get the latent dimension of the best ARI.
    best_latent_dim_ari = []
    for dataset in ari_res["Dataset"].unique():
        for method in ari_res["Method"].unique():
            for resolution in ari_res["Resolution"].unique()[:-1]:
                idx = (
                    (ari_res["Dataset"] == dataset)
                    & (ari_res["Method"] == method)
                    & (ari_res["Resolution"] == resolution)
                )
                if method == "Mowgli":
                    idx &= ari_res["H_rna regularization"] == 0.01
                    idx &= ari_res["H_atac regularization"] == 0.1
                    idx &= ari_res["H_adt regularization"] == 0.01
                    idx &= ari_res["w_regularization"] == 0.001
                best_id = ari_res.loc[idx, "max_ARI"].astype("float").argmax()
                best_latent_dim_ari.append(
                    {
                        "Dataset": dataset,
                        "Method": method,
                        "Latent dimension": ari_res.loc[idx, "Latent dimension"].iloc[
                            best_id
                        ],
                        "ARI": ari_res.loc[idx, "ARI"].iloc[best_id],
                        "Resolution": ari_res.loc[idx, "Resolution"].iloc[best_id],
                    }
                )
    best_latent_dim_ari = pd.DataFrame(best_latent_dim_ari)

    # Get the maximum purity score for each experiment.
    # Get the maximum Purity for each experiment.
    xp_mean_purity = purity_res.groupby("xp_name")["Purity score"].mean()

    # Add this as a column to the dataframe.
    purity_res["mean_Purity"] = purity_res["xp_name"].map(xp_mean_purity)

    best_latent_dim_purity = []
    for dataset in purity_res["Dataset"].unique():
        for method in purity_res["Method"].unique():
            for resolution in purity_res["k"].unique():
                idx = (
                    (purity_res["Dataset"] == dataset)
                    & (purity_res["Method"] == method)
                    & (purity_res["k"] == resolution)
                )
                if method == "Mowgli":
                    idx &= purity_res["H_rna regularization"] == 0.01
                    idx &= purity_res["H_atac regularization"] == 0.1
                    idx &= purity_res["H_adt regularization"] == 0.01
                    idx &= purity_res["w_regularization"] == 0.001
                best_id = purity_res.loc[idx, "mean_Purity"].astype("float").argmax()
                best_latent_dim_purity.append(
                    {
                        "Dataset": dataset,
                        "Method": method,
                        "Latent dimension": purity_res.loc[
                            idx, "Latent dimension"
                        ].iloc[best_id],
                        "Purity score": purity_res.loc[idx, "Purity score"].iloc[
                            best_id
                        ],
                        "k": purity_res.loc[idx, "k"].iloc[best_id],
                    }
                )
    best_latent_dim_purity = pd.DataFrame(best_latent_dim_purity)

    def plot_scores_hue_simulated(
        palette,
        method,
        path,
        hue="Latent dimension",
        hue_order=["5", "15", "30", "50"],
    ):
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

        # If the method if Mowgli, filter the regularization parameters.
        if method == "Mowgli":
            if hue != "Latent dimension":
                idx = idx & (scores_df["Latent dimension"] == "5")
            if hue != "H_rna regularization":
                idx = idx & (scores_df["H_rna regularization"] == 0.01)
            if hue != "H_atac regularization":
                idx = idx & (scores_df["H_atac regularization"] == 0.1)
            if hue != "H_adt regularization":
                idx = idx & (scores_df["H_adt regularization"] == 0.01)
            if hue != "w_regularization":
                idx = idx & (scores_df["w_regularization"] == 0.001)

        print(scores_df.loc[idx].index)

        sns.barplot(
            data=scores_df.loc[idx],
            y="Dataset",
            x="Silhouette score",
            hue=hue,
            hue_order=hue_order,
            palette=palette,
            order=datasets,
            ax=axes["A"],
        )
        axes["A"].legend(
            title=hue,
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

            # If the method if Mowgli, filter the regularization parameters.
            if method == "Mowgli":
                if hue != "Latent dimension":
                    idx = idx & (ari_res["Latent dimension"] == "5")
                if hue != "H_rna regularization":
                    idx = idx & (ari_res["H_rna regularization"] == 0.01)
                if hue != "H_atac regularization":
                    idx = idx & (ari_res["H_atac regularization"] == 0.1)
                if hue != "H_adt regularization":
                    idx = idx & (ari_res["H_adt regularization"] == 0.01)
                if hue != "w_regularization":
                    idx = idx & (ari_res["w_regularization"] == 0.001)

            print(datasets[i], idx.sum())
            sns.lineplot(
                data=ari_res.loc[idx],
                x="Resolution",
                y="ARI",
                hue=hue,
                hue_order=hue_order,
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

            # If the method if Mowgli, filter the regularization parameters.
            if method == "Mowgli":
                if hue != "Latent dimension":
                    idx = idx & (purity_res["Latent dimension"] == "5")
                if hue != "H_rna regularization":
                    idx = idx & (purity_res["H_rna regularization"] == 0.01)
                if hue != "H_atac regularization":
                    idx = idx & (purity_res["H_atac regularization"] == 0.1)
                if hue != "H_adt regularization":
                    idx = idx & (purity_res["H_adt regularization"] == 0.01)
                if hue != "w_regularization":
                    idx = idx & (purity_res["w_regularization"] == 0.001)

            sns.lineplot(
                data=purity_res.loc[idx],
                x="k",
                y="Purity score",
                hue=hue,
                hue_order=hue_order,
                palette=palette,
                ax=ax,
            )
            ax.get_legend().remove()
            ax.set(ylim=(0.33, ymax))
            if i < len(datasets) - 1:
                ax.set_xticklabels([])
                ax.set(xlabel=None)

        # Show a pretty grid.
        for i in axes:
            axes[i].grid()

        # Save the figure.
        fig.savefig(path)

    def plot_scores_hue_real(
        palette,
        method,
        path,
        hue="Latent dimension",
        hue_order=["5", "15", "30", "50"],
    ):

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

        # If the method is Mowgli, filter the regularization parameters.
        if method == "Mowgli":
            if hue != "Latent dimension":
                idx_liu = idx & scores_df["Dataset"].str.contains("Liu")
                idx_liu = idx_liu & (scores_df["Latent dimension"] == "5")

                idx_notliu = idx & ~scores_df["Dataset"].str.contains("Liu")
                idx_notliu = idx_notliu & (scores_df["Latent dimension"] == "50")

                idx = idx_liu | idx_notliu
            if hue != "H_rna regularization":
                idx = idx & (scores_df["H_rna regularization"] == 0.01)
            if hue != "H_atac regularization":
                idx = idx & (scores_df["H_atac regularization"] == 0.1)
            if hue != "H_adt regularization":
                idx = idx & (scores_df["H_adt regularization"] == 0.01)
            if hue != "w_regularization":
                idx = idx & (scores_df["w_regularization"] == 0.001)

        sns.barplot(
            data=scores_df.loc[idx],
            y="Dataset",
            x="Silhouette score",
            hue=hue,
            hue_order=hue_order,
            palette=palette,
            order=datasets,
            ax=axes["A"],
        )
        axes["A"].legend(
            title=hue,
            bbox_to_anchor=(-0.1, 1),
            ncol=1,
        )
        axes["A"].set_yticklabels(axes["A"].get_yticklabels(), rotation=90, va="center")

        # Visualize the ARI as a line plot.
        ymin = ari_res["ARI"].min()
        ymax = ari_res["ARI"].max() + 0.05
        for i, ax in enumerate([axes["B"], axes["D"], axes["F"], axes["H"], axes["J"]]):
            idx = ari_res["Dataset"].str.contains(datasets[i])
            idx = idx & (ari_res["Method"] == method)

            # If the method if Mowgli, filter the regularization parameters.
            if method == "Mowgli":
                if hue != "Latent dimension":
                    idx_liu = idx & ari_res["Dataset"].str.contains("Liu")
                    idx_liu = idx_liu & (ari_res["Latent dimension"] == "5")

                    idx_notliu = idx & ~ari_res["Dataset"].str.contains("Liu")
                    idx_notliu = idx_notliu & (ari_res["Latent dimension"] == "50")

                    idx = idx_liu | idx_notliu
                if hue != "H_rna regularization":
                    idx = idx & (ari_res["H_rna regularization"] == 0.01)
                if hue != "H_atac regularization":
                    idx = idx & (ari_res["H_atac regularization"] == 0.1)
                if hue != "H_adt regularization":
                    idx = idx & (ari_res["H_adt regularization"] == 0.01)
                if hue != "w_regularization":
                    idx = idx & (ari_res["w_regularization"] == 0.001)

            print(datasets[i], idx.sum())
            sns.lineplot(
                data=ari_res.loc[idx],
                x="Resolution",
                y="ARI",
                hue=hue,
                hue_order=hue_order,
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
        for i, ax in enumerate([axes["C"], axes["E"], axes["G"], axes["I"], axes["K"]]):
            idx = purity_res["Dataset"].str.contains(datasets[i])
            idx = idx & (purity_res["Method"] == method)

            # If the method if Mowgli, filter the regularization parameters.
            if method == "Mowgli":
                if hue != "Latent dimension":
                    idx_liu = idx & purity_res["Dataset"].str.contains("Liu")
                    idx_liu = idx_liu & (purity_res["Latent dimension"] == "5")

                    idx_notliu = idx & ~purity_res["Dataset"].str.contains("Liu")
                    idx_notliu = idx_notliu & (purity_res["Latent dimension"] == "50")

                    idx = idx_liu | idx_notliu
                if hue != "H_rna regularization":
                    idx = idx & (purity_res["H_rna regularization"] == 0.01)
                if hue != "H_atac regularization":
                    idx = idx & (purity_res["H_atac regularization"] == 0.1)
                if hue != "H_adt regularization":
                    idx = idx & (purity_res["H_adt regularization"] == 0.01)
                if hue != "w_regularization":
                    idx = idx & (purity_res["w_regularization"] == 0.001)

            sns.lineplot(
                data=purity_res.loc[idx],
                x="k",
                y="Purity score",
                hue=hue,
                hue_order=hue_order,
                palette=palette,
                ax=ax,
            )
            ax.get_legend().remove()
            ax.set(ylim=(0.60, ymax))
            if i < len(datasets) - 1:
                ax.set_xticklabels([])
                ax.set(xlabel=None)

        # Show a pretty grid.
        for i in axes:
            axes[i].grid()

        # Save the figure.
        fig.savefig(path)

    def plot_best_scores_real(path):

        dataset_labels = [
            "Liu",
            "10X PBMC",
            "Open Problems Multiome",
            "Open Problems CITE-seq",
            "Bone Marrow CITE-seq",
        ]
        fig = plt.figure(constrained_layout=True, figsize=(15, 20))
        axes = fig.subplot_mosaic(
            """
            ABC
            ADE
            AFG
            AHI
            AJK
            """
        )

        # Plot silhouette.
        method_labels = ["Seurat", "MOFA+", "NMF", "Cobolt", "Multigrate", "Mowgli"]
        sns.barplot(
            data=best_latent_dim_sil,
            y="Dataset",
            x="Silhouette score",
            hue="Method",
            hue_order=method_labels,
            order=dataset_labels,
            ax=axes["A"],
        )
        axes["A"].set_yticklabels(
            axes["A"].get_yticklabels(), rotation=90, verticalalignment="center"
        )
        for i, container in enumerate(axes["A"].containers):
            dim_labels = []
            for label in dataset_labels:
                idx = best_latent_dim_sil["Method"] == method_labels[i]
                idx &= best_latent_dim_sil["Dataset"] == label
                dim_labels.append(
                    best_latent_dim_sil.loc[idx, "Latent dimension"].values[0]
                )
            axes["A"].bar_label(container, labels=dim_labels)

        # Plot ARI.
        ari_axes = ["B", "D", "F", "H", "J", "L"]
        ymin = 0.25#best_latent_dim_ari["ARI"].min()
        ymax = best_latent_dim_ari["ARI"].max() + 0.05
        for i_dataset, dataset in enumerate(dataset_labels):
            sns.lineplot(
                data=best_latent_dim_ari[best_latent_dim_ari["Dataset"] == dataset],
                x="Resolution",
                y="ARI",
                hue="Method",
                hue_order=method_labels,
                legend=None,
                ax=axes[ari_axes[i_dataset]],
            )
            axes[ari_axes[i_dataset]].set(ylim=(ymin, ymax))
            # Add the latent dimension in the legend.
            for i, line in enumerate(axes[ari_axes[i_dataset]].lines):
                latent_dim = best_latent_dim_ari.loc[
                    (best_latent_dim_ari["Dataset"] == dataset)
                    & (best_latent_dim_ari["Method"] == method_labels[i]),
                    "Latent dimension",
                ].values[0]
                if method_labels[i] == "Seurat":
                    line.set_label(f"{method_labels[i]}")
                else:
                    line.set_label(f"{method_labels[i]} (latent dimensions: {latent_dim})")
            axes[ari_axes[i_dataset]].legend()

        # Plot purity.
        purity_axes = ["C", "E", "G", "I", "K", "M"]
        ymin = 0.65#best_latent_dim_purity["Purity score"].min()
        ymax = 1
        for i_dataset, dataset in enumerate(dataset_labels):
            sns.lineplot(
                data=best_latent_dim_purity[
                    best_latent_dim_purity["Dataset"] == dataset
                ],
                x="k",
                y="Purity score",
                hue="Method",
                hue_order=method_labels,
                legend=None,
                ax=axes[purity_axes[i_dataset]],
            )
            axes[purity_axes[i_dataset]].set(ylim=(ymin, ymax))
            # Add the latent dimension in the legend.
            for i, line in enumerate(axes[purity_axes[i_dataset]].lines):
                latent_dim = best_latent_dim_purity.loc[
                    (best_latent_dim_purity["Dataset"] == dataset)
                    & (best_latent_dim_purity["Method"] == method_labels[i]),
                    "Latent dimension",
                ].values[0]
                if method_labels[i] == "Seurat":
                    line.set_label(f"{method_labels[i]}")
                else:
                    line.set_label(f"{method_labels[i]} (latent dimensions: {latent_dim})")
            axes[purity_axes[i_dataset]].legend()

        for i in axes:
            axes[i].grid()

        plt.savefig(path)

    def plot_best_scores_simulated(path):
        dataset_labels = [
            "Simulated 1",
            "Simulated 2",
            "Simulated 3",
            "Simulated 4",
            "Simulated 5",
            "Simulated 6",
        ]
        fig = plt.figure(constrained_layout=True, figsize=(15, 20))
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

        # Plot silhouette.
        method_labels = ["Seurat", "MOFA+", "NMF", "Cobolt", "Multigrate", "Mowgli"]
        sns.barplot(
            data=best_latent_dim_sil,
            y="Dataset",
            x="Silhouette score",
            hue="Method",
            hue_order=method_labels,
            order=dataset_labels,
            ax=axes["A"],
        )
        axes["A"].set_yticklabels(
            axes["A"].get_yticklabels(), rotation=90, verticalalignment="center"
        )
        for i, container in enumerate(axes["A"].containers):
            dim_labels = []
            for label in dataset_labels:
                idx = best_latent_dim_sil["Method"] == method_labels[i]
                idx &= best_latent_dim_sil["Dataset"] == label
                dim_labels.append(
                    best_latent_dim_sil.loc[idx, "Latent dimension"].values[0]
                )
            axes["A"].bar_label(container, labels=dim_labels)

        # Plot ARI.
        ari_axes = ["B", "D", "F", "H", "J", "L"]
        for i_dataset, dataset in enumerate(dataset_labels):
            sns.lineplot(
                data=best_latent_dim_ari[best_latent_dim_ari["Dataset"] == dataset],
                x="Resolution",
                y="ARI",
                hue="Method",
                hue_order=method_labels,
                legend=None,
                ax=axes[ari_axes[i_dataset]],
            )
            # Add the latent dimension in the legend.
            for i, line in enumerate(axes[ari_axes[i_dataset]].lines):
                latent_dim = best_latent_dim_ari.loc[
                    (best_latent_dim_ari["Dataset"] == dataset)
                    & (best_latent_dim_ari["Method"] == method_labels[i]),
                    "Latent dimension",
                ].values[0]
                if method_labels[i] == "Seurat":
                    line.set_label(f"{method_labels[i]}")
                else:
                    line.set_label(f"{method_labels[i]} (latent dimensions: {latent_dim})")
            axes[ari_axes[i_dataset]].legend()

        # Plot purity.
        purity_axes = ["C", "E", "G", "I", "K", "M"]
        for i_dataset, dataset in enumerate(dataset_labels):
            sns.lineplot(
                data=best_latent_dim_purity[
                    best_latent_dim_purity["Dataset"] == dataset
                ],
                x="k",
                y="Purity score",
                hue="Method",
                hue_order=method_labels,
                legend=None,
                ax=axes[purity_axes[i_dataset]],
            )
            # Add the latent dimension in the legend.
            for i, line in enumerate(axes[purity_axes[i_dataset]].lines):
                latent_dim = best_latent_dim_purity.loc[
                    (best_latent_dim_purity["Dataset"] == dataset)
                    & (best_latent_dim_purity["Method"] == method_labels[i]),
                    "Latent dimension",
                ].values[0]
                if method_labels[i] == "Seurat":
                    line.set_label(f"{method_labels[i]}")
                else:
                    line.set_label(f"{method_labels[i]} (latent dimensions: {latent_dim})")
            axes[purity_axes[i_dataset]].legend()

        for i in axes:
            axes[i].grid()

        plt.savefig(path)

    def plot_scores_selected_simulated(path):

        # Select only certain latent dimensions
        idx = (
            (scores_df["Latent dimension"] == "5")
            & (scores_df["Method"] == "MOFA+")
            & (scores_df["Dataset"].str.contains("Simulated"))
        )
        idx |= (
            (scores_df["Latent dimension"] == "5")
            & (scores_df["Method"] == "Mowgli")
            & (scores_df["H_rna regularization"] == 0.01)
            & (scores_df["H_atac regularization"] == 0.1)
            & (scores_df["w_regularization"] == 0.001)
            & (scores_df["Dataset"].str.contains("Simulated"))
        )
        idx |= (
            (scores_df["Latent dimension"] == "5")
            & (scores_df["Method"] == "NMF")
            & (scores_df["Dataset"].str.contains("Simulated"))
        )
        idx |= (
            (scores_df["Latent dimension"] == "5")
            & (scores_df["Method"] == "Cobolt")
            & (scores_df["Dataset"].str.contains("Simulated"))
        )
        idx |= (
            (scores_df["Latent dimension"] == "15")
            & (scores_df["Method"] == "Multigrate")
            & (scores_df["Dataset"].str.contains("Simulated"))
        )
        idx |= scores_df["Method"] == "Seurat"
        scores_df_sub = scores_df.loc[idx]

        idx = (
            (ari_res["Latent dimension"] == "5")
            & (ari_res["Method"] == "MOFA+")
            & (ari_res["Dataset"].str.contains("Simulated"))
        )
        idx |= (
            (ari_res["Latent dimension"] == "5")
            & (ari_res["Method"] == "Mowgli")
            & (ari_res["H_rna regularization"] == 0.01)
            & (ari_res["H_atac regularization"] == 0.1)
            & (ari_res["w_regularization"] == 0.001)
            & (ari_res["Dataset"].str.contains("Simulated"))
        )
        idx |= (
            (ari_res["Latent dimension"] == "5")
            & (ari_res["Method"] == "NMF")
            & (ari_res["Dataset"].str.contains("Simulated"))
        )
        idx |= (
            (ari_res["Latent dimension"] == "5")
            & (ari_res["Method"] == "Cobolt")
            & (ari_res["Dataset"].str.contains("Simulated"))
        )
        idx |= (
            (ari_res["Latent dimension"] == "15")
            & (ari_res["Method"] == "Multigrate")
            & (ari_res["Dataset"].str.contains("Simulated"))
        )
        idx |= ari_res["Method"] == "Seurat"
        ari_res_sub = ari_res.loc[idx]

        idx = (
            (purity_res["Latent dimension"] == "5")
            & (purity_res["Method"] == "MOFA+")
            & (purity_res["Dataset"].str.contains("Simulated"))
        )
        idx |= (
            (purity_res["Latent dimension"] == "5")
            & (purity_res["Method"] == "Mowgli")
            & (purity_res["H_rna regularization"] == 0.01)
            & (purity_res["H_atac regularization"] == 0.1)
            & (purity_res["w_regularization"] == 0.001)
            & (purity_res["Dataset"].str.contains("Simulated"))
        )
        idx |= (
            (purity_res["Latent dimension"] == "5")
            & (purity_res["Method"] == "NMF")
            & (purity_res["Dataset"].str.contains("Simulated"))
        )
        idx |= (
            (purity_res["Latent dimension"] == "5")
            & (purity_res["Method"] == "Cobolt")
            & (purity_res["Dataset"].str.contains("Simulated"))
        )
        idx |= (
            (purity_res["Latent dimension"] == "15")
            & (purity_res["Method"] == "Multigrate")
            & (purity_res["Dataset"].str.contains("Simulated"))
        )
        idx |= purity_res["Method"] == "Seurat"
        purity_res_sub = purity_res.loc[idx]

        # List the datasets.
        datasets = [f"Simulated {i}" for i in range(1, 7)]

        # Define the order of the datasets.
        for i in range(len(datasets)):
            idx = scores_df_sub["Dataset"].str.contains(datasets[i])
            scores_df_sub.loc[idx, "Order"] = i

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
        idx = scores_df_sub["Dataset"].str.contains("|".join(datasets))
        sns.barplot(
            data=scores_df_sub.loc[idx],
            y="Dataset",
            x="Silhouette score",
            hue="Method",
            hue_order=["Seurat", "MOFA+", "NMF", "Cobolt", "Multigrate", "Mowgli"],
            order=datasets,
            ax=axes["A"],
        )
        axes["A"].legend(
            title="Method",
            bbox_to_anchor=(-0.1, 1),
            ncol=1,
        )
        axes["A"].set_yticklabels(
            [
                "Mix groups in RNA",
                "Mix in both omics",
                "Rare population",
                "82\% sparse",
                "90\% sparse",
                "96\% sparse",
            ],
            rotation=90,
            va="center",
        )

        # Visualize the ARI as a line plot.
        ymin = ari_res_sub["ARI"].min()
        ymax = ari_res_sub["ARI"].max() + 0.05
        for i, ax in enumerate(
            [axes["B"], axes["D"], axes["F"], axes["H"], axes["J"], axes["L"]]
        ):
            idx = ari_res_sub["Dataset"].str.contains(datasets[i])
            sns.barplot(
                data=ari_res_sub.loc[idx],
                x="Dataset",
                y="ARI",
                hue="Method",
                hue_order=["Seurat", "MOFA+", "NMF", "Cobolt", "Multigrate", "Mowgli"],
                estimator=max,
                ci=None,
                ax=ax,
            )
            ax.get_legend().remove()
            ax.set(ylim=(ymin, ymax))
            ax.set(ylabel="maximum ARI")
            ax.set_xticks([])
            if i < len(datasets) - 1:
                ax.set(xlabel=None)

        # Visualize the purity score as a line plot.
        ymin = purity_res_sub["Purity score"].min()
        ymax = purity_res_sub["Purity score"].max() + 0.05
        for i, ax in enumerate(
            [axes["C"], axes["E"], axes["G"], axes["I"], axes["K"], axes["M"]]
        ):
            idx = purity_res_sub["Dataset"].str.contains(datasets[i])
            sns.lineplot(
                data=purity_res_sub.loc[idx],
                x="k",
                y="Purity score",
                hue="Method",
                hue_order=["Seurat", "MOFA+", "NMF", "Cobolt", "Multigrate", "Mowgli"],
                ax=ax,
            )
            ax.get_legend().remove()
            ax.set(ylim=(0.45, ymax))
            ax.set(xlim=(0, 20))
            if i < len(datasets) - 1:
                ax.set_xticklabels([])
                ax.set(xlabel=None)

        for i in axes:
            axes[i].grid()

        plt.savefig(path)

    def plot_scores_selected_real(path):

        # Select only certain latent dimensions
        idx = (
            (scores_df["Latent dimension"] == "5")
            & (scores_df["Method"] == "MOFA+")
            & (scores_df["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (scores_df["Latent dimension"] == "15")
            & (scores_df["Method"] == "MOFA+")
            & (~scores_df["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (scores_df["Latent dimension"] == "5")
            & (scores_df["Method"] == "Mowgli")
            & (scores_df["H_rna regularization"] == 0.01)
            & (scores_df["H_atac regularization"] == 0.1)
            & (scores_df["H_adt regularization"] == 0.01)
            & (scores_df["w_regularization"] == 0.001)
            & (scores_df["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (scores_df["Latent dimension"] == "50")
            & (scores_df["Method"] == "Mowgli")
            & (scores_df["H_rna regularization"] == 0.01)
            & (scores_df["H_atac regularization"] == 0.1)
            & (scores_df["H_adt regularization"] == 0.01)
            & (scores_df["w_regularization"] == 0.001)
            & (~scores_df["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (scores_df["Latent dimension"] == "5")
            & (scores_df["Method"] == "NMF")
            & (scores_df["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (scores_df["Latent dimension"] == "30")
            & (scores_df["Method"] == "NMF")
            & (~scores_df["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (scores_df["Latent dimension"] == "5")
            & (scores_df["Method"] == "Cobolt")
            & (scores_df["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (scores_df["Latent dimension"] == "30")
            & (scores_df["Method"] == "Cobolt")
            & (~scores_df["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (scores_df["Latent dimension"] == "15")
            & (scores_df["Method"] == "Multigrate")
            & (scores_df["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (scores_df["Latent dimension"] == "50")
            & (scores_df["Method"] == "Multigrate")
            & (~scores_df["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= scores_df["Method"] == "Seurat"
        scores_df_sub = scores_df.loc[idx]

        idx = (
            (ari_res["Latent dimension"] == "5")
            & (ari_res["Method"] == "MOFA+")
            & (ari_res["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (ari_res["Latent dimension"] == "15")
            & (ari_res["Method"] == "MOFA+")
            & (~ari_res["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (ari_res["Latent dimension"] == "5")
            & (ari_res["Method"] == "Mowgli")
            & (ari_res["H_rna regularization"] == 0.01)
            & (ari_res["H_atac regularization"] == 0.1)
            & (ari_res["H_adt regularization"] == 0.01)
            & (ari_res["w_regularization"] == 0.001)
            & (ari_res["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (ari_res["Latent dimension"] == "50")
            & (ari_res["Method"] == "Mowgli")
            & (ari_res["H_rna regularization"] == 0.01)
            & (ari_res["H_atac regularization"] == 0.1)
            & (ari_res["H_adt regularization"] == 0.01)
            & (ari_res["w_regularization"] == 0.001)
            & (~ari_res["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (ari_res["Latent dimension"] == "5")
            & (ari_res["Method"] == "NMF")
            & (ari_res["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (ari_res["Latent dimension"] == "30")
            & (ari_res["Method"] == "NMF")
            & (~ari_res["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (ari_res["Latent dimension"] == "5")
            & (ari_res["Method"] == "Cobolt")
            & (ari_res["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (ari_res["Latent dimension"] == "30")
            & (ari_res["Method"] == "Cobolt")
            & (~ari_res["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (ari_res["Latent dimension"] == "15")
            & (ari_res["Method"] == "Multigrate")
            & (ari_res["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (ari_res["Latent dimension"] == "50")
            & (ari_res["Method"] == "Multigrate")
            & (~ari_res["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= ari_res["Method"] == "Seurat"
        ari_res_sub = ari_res.loc[idx]

        idx = (
            (purity_res["Latent dimension"] == "5")
            & (purity_res["Method"] == "MOFA+")
            & (purity_res["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (purity_res["Latent dimension"] == "15")
            & (purity_res["Method"] == "MOFA+")
            & (~purity_res["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (purity_res["Latent dimension"] == "5")
            & (purity_res["Method"] == "Mowgli")
            & (purity_res["H_rna regularization"] == 0.01)
            & (purity_res["H_atac regularization"] == 0.1)
            & (purity_res["H_adt regularization"] == 0.01)
            & (purity_res["w_regularization"] == 0.001)
            & (purity_res["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (purity_res["Latent dimension"] == "50")
            & (purity_res["Method"] == "Mowgli")
            & (purity_res["H_rna regularization"] == 0.01)
            & (purity_res["H_atac regularization"] == 0.1)
            & (purity_res["H_adt regularization"] == 0.01)
            & (purity_res["w_regularization"] == 0.001)
            & (~purity_res["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (purity_res["Latent dimension"] == "5")
            & (purity_res["Method"] == "NMF")
            & (purity_res["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (purity_res["Latent dimension"] == "30")
            & (purity_res["Method"] == "NMF")
            & (~purity_res["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (purity_res["Latent dimension"] == "5")
            & (purity_res["Method"] == "Cobolt")
            & (purity_res["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (purity_res["Latent dimension"] == "30")
            & (purity_res["Method"] == "Cobolt")
            & (~purity_res["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (purity_res["Latent dimension"] == "15")
            & (purity_res["Method"] == "Multigrate")
            & (purity_res["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= (
            (purity_res["Latent dimension"] == "50")
            & (purity_res["Method"] == "Multigrate")
            & (~purity_res["Dataset"].str.contains("Liu|Simulated"))
        )
        idx |= purity_res["Method"] == "Seurat"
        purity_res_sub = purity_res.loc[idx]

        # List the datasets.
        datasets = [
            "Liu",
            "10X PBMC",
            "Open Problems Multiome",
            "Open Problems CITE-seq",
            "Bone Marrow CITE-seq",
        ]

        for i in range(len(datasets)):
            idx = scores_df_sub["Dataset"].str.contains(datasets[i])
            scores_df_sub.loc[idx, "Order"] = i

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
        idx = scores_df_sub["Dataset"].str.contains("|".join(datasets))
        sns.barplot(
            data=scores_df_sub.loc[idx],
            y="Dataset",
            x="Silhouette score",
            hue="Method",
            hue_order=["Seurat", "MOFA+", "NMF", "Cobolt", "Multigrate", "Mowgli"],
            order=datasets,
            ax=axes["A"],
        )
        axes["A"].legend(
            title="Method",
            bbox_to_anchor=(-0.1, 1),
            ncol=1,
        )
        axes["A"].set_yticklabels(axes["A"].get_yticklabels(), rotation=90, va="center")

        # Visualize the ARI as a line plot.
        ymin = ari_res_sub["ARI"].min()
        ymax = ari_res_sub["ARI"].max() + 0.05
        for i, ax in enumerate([axes["B"], axes["D"], axes["F"], axes["H"], axes["J"]]):
            idx = ari_res_sub["Dataset"].str.contains(datasets[i])
            sns.lineplot(
                data=ari_res_sub.loc[idx],
                x="Resolution",
                y="ARI",
                hue="Method",
                hue_order=["Seurat", "MOFA+", "NMF", "Cobolt", "Multigrate", "Mowgli"],
                ax=ax,
            )
            ax.get_legend().remove()
            ax.set(ylim=(ymin, ymax))
            if i < len(datasets) - 1:
                ax.set_xticklabels([])
                ax.set(xlabel=None)

        # Visualize the purity score as a line plot.
        ymin = purity_res_sub["Purity score"].min()
        ymax = purity_res_sub["Purity score"].max() + 0.05
        for i, ax in enumerate([axes["C"], axes["E"], axes["G"], axes["I"], axes["K"]]):
            idx = purity_res_sub["Dataset"].str.contains(datasets[i])
            sns.lineplot(
                data=purity_res_sub.loc[idx],
                x="k",
                y="Purity score",
                hue="Method",
                hue_order=["Seurat", "MOFA+", "NMF", "Cobolt", "Multigrate", "Mowgli"],
                ax=ax,
            )
            ax.get_legend().remove()
            ax.set(xlim=(0, 20))
            ax.set(ylim=(0.7, ymax))
            if i < len(datasets) - 1:
                ax.set_xticklabels([])
                ax.set(xlabel=None)

        for i in axes:
            axes[i].grid()

        plt.savefig(path)

    # Plot the scores for Mowgli for simulated data, per dimension.
    print("Plotting scores for Mowgli for simulated data, per dimension.")
    path = f"{cfg.figure_path}scores_dim_simulated_mowgli.pdf"
    palette = sns.color_palette("flare", n_colors=5)
    plot_scores_hue_simulated(
        palette,
        "Mowgli",
        path,
        hue="Latent dimension",
        hue_order=["5", "15", "30", "50", "100"],
    )

    # Plot the scores for Mowgli for simulated data, per H_rna regularization.
    print("Plotting scores for Mowgli for simulated data, per H_rna regularization.")
    path = f"{cfg.figure_path}scores_h_rna_simulated_mowgli.pdf"
    palette = sns.color_palette("flare", n_colors=4)
    plot_scores_hue_simulated(
        palette,
        "Mowgli",
        path,
        hue="H_rna regularization",
        hue_order=[1.0, 0.1, 0.01, 0.001],
    )

    # Plot the scores for Mowgli for simulated data, per H_atac regularization.
    print("Plotting scores for Mowgli for simulated data, per H_atac regularization.")
    path = f"{cfg.figure_path}scores_h_atac_simulated_mowgli.pdf"
    palette = sns.color_palette("flare", n_colors=4)
    plot_scores_hue_simulated(
        palette,
        "Mowgli",
        path,
        hue="H_atac regularization",
        hue_order=[1.0, 0.1, 0.01, 0.001],
    )

    # Plot the scores for Mowgli for simulated data, per w_regularization.
    print("Plotting scores for Mowgli for simulated data, per w_regularization.")
    path = f"{cfg.figure_path}scores_w_regularization_simulated_mowgli.pdf"
    palette = sns.color_palette("flare", n_colors=3)
    plot_scores_hue_simulated(
        palette,
        "Mowgli",
        path,
        hue="w_regularization",
        hue_order=[0.01, 0.001, 0.0001],
    )

    # Plot the scores for cobolt for simulated data, per dimension.
    print("Plotting scores for cobolt for simulated data, per dimension.")
    path = f"{cfg.figure_path}scores_dim_simulated_cobolt.pdf"
    palette = sns.color_palette("flare", n_colors=5)
    plot_scores_hue_simulated(
        palette,
        "Cobolt",
        path,
        hue="Latent dimension",
        hue_order=["5", "15", "30", "50", "100"],
    )

    # Plot the scores for multigrate for simulated data, per dimension.
    print("Plotting scores for multigrate for simulated data, per dimension.")
    path = f"{cfg.figure_path}scores_dim_simulated_multigrate.pdf"
    palette = sns.color_palette("flare", n_colors=5)
    plot_scores_hue_simulated(
        palette,
        "Multigrate",
        path,
        hue="Latent dimension",
        hue_order=["5", "15", "30", "50", "100"],
    )

    # Plot the scores for MOFA+ for simulated data, per dimension.
    print("Plotting scores for MOFA+ for simulated data, per dimension.")
    path = f"{cfg.figure_path}scores_dim_simulated_mofa.pdf"
    palette = sns.color_palette("flare", n_colors=4)
    plot_scores_hue_simulated(palette, "MOFA+", path, hue="Latent dimension")

    # Plot the scores for NMF for simulated data, per dimension.
    print("Plotting scores for NMF for simulated data, per dimension.")
    path = f"{cfg.figure_path}scores_dim_simulated_nmf.pdf"
    palette = sns.color_palette("flare", n_colors=4)
    plot_scores_hue_simulated(palette, "NMF", path, hue="Latent dimension")

    # Plot the scores for Mowgli for real data, per dimension.
    print("Plotting scores for Mowgli for real data, per dimension.")
    path = f"{cfg.figure_path}scores_dim_real_mowgli.pdf"
    palette = sns.color_palette("flare", n_colors=5)
    plot_scores_hue_real(
        palette, "Mowgli", path, hue_order=["5", "15", "30", "50", "100"]
    )

    # Plot the scores for Mowgli for real data, per H_rna regularization.
    print("Plotting scores for Mowgli for real data, per H_rna regularization.")
    path = f"{cfg.figure_path}scores_h_rna_real_mowgli.pdf"
    palette = sns.color_palette("flare", n_colors=4)
    plot_scores_hue_real(
        palette,
        "Mowgli",
        path,
        hue="H_rna regularization",
        hue_order=[1.0, 0.1, 0.01, 0.001],
    )

    # Plot the scores for Mowgli for real data, per H_atac regularization.
    print("Plotting scores for Mowgli for real data, per H_atac regularization.")
    path = f"{cfg.figure_path}scores_h_atac_real_mowgli.pdf"
    palette = sns.color_palette("flare", n_colors=4)
    plot_scores_hue_real(
        palette,
        "Mowgli",
        path,
        hue="H_atac regularization",
        hue_order=[1.0, 0.1, 0.01, 0.001],
    )

    # Plot the scores for Mowgli for real data, per H_adt regularization.
    print("Plotting scores for Mowgli for real data, per H_adt regularization.")
    path = f"{cfg.figure_path}scores_h_adt_real_mowgli.pdf"
    palette = sns.color_palette("flare", n_colors=4)
    plot_scores_hue_real(
        palette,
        "Mowgli",
        path,
        hue="H_adt regularization",
        hue_order=[1.0, 0.1, 0.01, 0.001],
    )

    # Plot the scores for Mowgli for real data, per w_regularization.
    print("Plotting scores for Mowgli for real data, per w_regularization.")
    path = f"{cfg.figure_path}scores_w_regularization_real_mowgli.pdf"
    palette = sns.color_palette("flare", n_colors=3)
    plot_scores_hue_real(
        palette,
        "Mowgli",
        path,
        hue="w_regularization",
        hue_order=[0.01, 0.001, 0.0001],
    )

    # Plot the scores for cobolt for real data, per dimension.
    print("Plotting scores for Cobolt for real data, per dimension.")
    path = f"{cfg.figure_path}scores_dim_real_cobolt.pdf"
    palette = sns.color_palette("flare", n_colors=5)
    plot_scores_hue_real(
        palette,
        "Cobolt",
        path,
        hue="Latent dimension",
        hue_order=["5", "15", "30", "50", "100"],
    )

    # Plot the scores for multigrate for real data, per dimension.
    print("Plotting scores for multigrate for real data, per dimension.")
    path = f"{cfg.figure_path}scores_dim_real_multigrate.pdf"
    palette = sns.color_palette("flare", n_colors=5)
    plot_scores_hue_real(
        palette,
        "Multigrate",
        path,
        hue="Latent dimension",
        hue_order=["5", "15", "30", "50", "100"],
    )

    # Plot the scores for MOFA+ for real data, per dimension.
    print("Plotting scores for MOFA+ for real data, per dimension.")
    path = f"{cfg.figure_path}scores_dim_real_mofa.pdf"
    palette = sns.color_palette("flare", n_colors=4)
    plot_scores_hue_real(palette, "MOFA+", path, hue="Latent dimension")

    # Plot the scores for NMF for real data, per dimension.
    print("Plotting scores for NMF for real data, per dimension.")
    path = f"{cfg.figure_path}scores_dim_real_nmf.pdf"
    palette = sns.color_palette("flare", n_colors=4)
    plot_scores_hue_real(palette, "NMF", path, hue="Latent dimension")

    # Plot the scores for all methods for simulated data.
    print("Plotting scores for all methods for simulated data.")
    path = f"{cfg.figure_path}scores_simulated.pdf"
    palette = tue_cycler.cycler(color=tue_palettes.high_contrast)
    plot_scores_selected_simulated(path)

    # Plot the scores for all methods for real data.
    print("Plotting scores for all methods for real data.")
    path = f"{cfg.figure_path}scores_real.pdf"
    palette = tue_cycler.cycler(color=tue_palettes.high_contrast)
    plot_scores_selected_real(path)

    path = f"{cfg.figure_path}best_scores_real.pdf"
    plot_best_scores_real(path)

    path = f"{cfg.figure_path}best_scores_simulated.pdf"
    plot_best_scores_simulated(path)


if __name__ == "__main__":
    my_app()
