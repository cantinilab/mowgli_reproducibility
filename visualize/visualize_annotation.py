import hydra
from omegaconf import DictConfig


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def my_app(cfg: DictConfig) -> None:

    ####################################### Imports ######################################

    import anndata as ad
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    import mowgli
    import muon as mu
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import seaborn as sns
    from sklearn.metrics import silhouette_samples

    plt.rcParams["font.family"] = "Roboto Condensed"
    sc.set_figure_params(vector_friendly=True, dpi_save=300)

    ####################################### BMCITE #######################################
    # Load the data
    mdata = mu.read_h5mu(cfg.data_path + "BMCITE/bmcite_preprocessed.h5mu")

    # Load the Mowgli embedding
    path = cfg.w_path + "bmcite_mowgli_cosine_50_0_05_0_01_0_001.npy"
    mdata.obsm["X_mowgli"] = np.load(path, allow_pickle=True).item()["W"]

    # Load the Mowgli loadings
    path = cfg.h_path + "bmcite_mowgli_cosine_50_0_05_0_01_0_001.npy"
    H_OT = np.load(path, allow_pickle=True).item()
    mdata["rna"].uns["H_OT"] = H_OT["H_rna"]
    mdata["adt"].uns["H_OT"] = H_OT["H_adt"]

    # Compute the UMAP emebedding
    sc.pp.neighbors(mdata, use_rep="X_mowgli", n_neighbors=30)
    sc.tl.umap(mdata)

    # Compute the silhouette scores
    mdata.obs["silhouette"] = silhouette_samples(
        mdata.obsm["X_mowgli"],
        mdata.obs["rna:celltype"],
    )

    # Plot the cell-wise silhouette score on the UMAP
    sc.pl.umap(mdata, color=["silhouette", "rna:celltype"], show=False)
    plt.savefig(cfg.figure_path + "bmcite_silhouette.pdf", bbox_inches="tight")

    # Focus on two subsets of CD8+ T cells
    mdata = mdata[mdata.obs["rna:celltype"].str.contains("CD8 Effector_1|CD8 Memory_2")]
    mdata.obs["rna:celltype"] = pd.Categorical(mdata.obs["rna:celltype"].tolist())

    # Clustermap parameters
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    handles = [
        mpatches.Patch(color=colors[i], label=c)
        for i, c in enumerate(mdata.obs["rna:celltype"].cat.categories)
    ]
    clustermap_kwds = {
        "cmap": "viridis",
        "obs_keys": "celltype",
        "figsize": (5, 5),
        "col_cluster": False,
        "colors_ratio": 0.1,
        "yticklabels": False,
        "dendrogram_ratio": 0.1,
        "show": False,
    }

    # Plot a clustermap of ADT signal
    mdata["adt"].obs["celltype"] = mdata.obs["rna:celltype"]
    g = sc.pl.clustermap(mdata["adt"], cbar_pos=(1, 0.2, 0.025, 0.6), **clustermap_kwds)
    g.ax_heatmap.set_xticks(np.arange(0.5, len(mdata["adt"].var_names)))
    g.ax_heatmap.set_xticklabels(mdata["adt"].var_names)
    g.ax_heatmap.legend(
        handles=handles,
        bbox_to_anchor=(1.1, 1.1),
        loc="upper left",
        borderaxespad=0.0,
        frameon=False,
    )
    g.ax_heatmap.set_title("Hierarchical clustering of ADT signal", fontsize=14)
    plt.savefig(cfg.figure_path + "bmcite_adt_clustermap.pdf", bbox_inches="tight")

    # Plot a clustermap of the joint embedding
    embed = ad.AnnData(mdata.obsm["X_mowgli"], obs=mdata.obs)
    embed.obs["celltype"] = mdata.obs["rna:celltype"]
    g = sc.pl.clustermap(embed, cbar_pos=(1, 0.1, 0.025, 0.7), **clustermap_kwds)
    g.ax_heatmap.set_xticks(np.arange(0.5, 50, 5))
    g.ax_heatmap.set_xticklabels(np.arange(0, 50, 5))
    g.ax_heatmap.set_xticks(np.arange(0.5, 50, 1), minor=True)
    g.ax_heatmap.legend(
        handles=handles,
        bbox_to_anchor=(1.1, 1.1),
        loc="upper left",
        borderaxespad=0.0,
        frameon=False,
    )
    g.ax_heatmap.set_title("Hierarchical clustering of Mowgli embeddings", fontsize=14)
    plt.savefig(cfg.figure_path + "bmcite_mowgli_clustermap.pdf", bbox_inches="tight")

    # Plot the compared factors
    mdata.obs[mdata["adt"].var_names] = mdata["adt"].X
    xx, yy = mdata.obs["adt:CD45RO"], mdata.obs["adt:CD57"]
    fig, axes = plt.subplots(2, 4, figsize=(20, 8), tight_layout=True)

    sns.scatterplot(x=xx, y=yy, hue=mdata.obs["rna:celltype"], ax=axes[0, 0], alpha=0.7)
    axes[0, 0].set_title("Protein expression for two CD8 subtypes")

    sc.pl.stacked_violin(
        mdata["adt"],
        [
            "adt:CD8a",
            "adt:CD3",
            "adt:CD69",
            "adt:CD45RO",
            "adt:CD57",
            "adt:CD45RA",
            "adt:CD127-IL7Ra",
            "adt:CD27",
            "adt:CD28",
            "adt:CD56",
        ],
        groupby="celltype",
        ax=axes[1, 0],
        use_raw=False,
        show=False,
    )

    hue = mdata.obsm["X_mowgli"][:, 44]
    sns.scatterplot(x=xx, y=yy, hue=hue, palette="viridis", ax=axes[0, 1], alpha=0.7)
    axes[0, 1].set_title("Weight of factor 44 across cells")

    plt.sca(axes[1, 1])
    mowgli.pl.top_features(mdata, dim=44, mod="adt")
    axes[1, 1].set_title("Top ADT weights of factor 44")

    hue = mdata.obsm["X_mowgli"][:, 22]
    sns.scatterplot(x=xx, y=yy, hue=hue, palette="viridis", ax=axes[0, 2], alpha=0.7)
    axes[0, 2].set_title("Weight of factor 22 across cells")

    plt.sca(axes[1, 2])
    mowgli.pl.top_features(mdata, dim=22, mod="adt")
    axes[1, 2].set_title("Top ADT weights of factor 22")

    hue = mdata.obsm["X_mowgli"][:, 40]
    sns.scatterplot(x=xx, y=yy, hue=hue, palette="viridis", ax=axes[0, 3], alpha=0.7)
    axes[0, 3].set_title("Weight of factor 40 across cells")

    plt.sca(axes[1, 3])
    mowgli.pl.top_features(mdata, dim=40, mod="adt")
    axes[1, 3].set_title("Top ADT weights of factor 40")

    plt.savefig(cfg.figure_path + "bmcite_recap.pdf", bbox_inches="tight")

    ####################################### OPCITE #######################################

    # Load the data
    mdata = mu.read_h5mu(cfg.data_path + "OPCITE/opcite_preprocessed.h5mu.gz")

    # Load the Mowgli embedding
    path = cfg.w_path + "opcite_mowgli_cosine_50_0_05_0_01_0_001.npy"
    mdata.obsm["X_mowgli"] = np.load(path, allow_pickle=True).item()["W"]

    # Load the Mowgli loadings
    path = cfg.h_path + "opcite_mowgli_cosine_50_0_05_0_01_0_001.npy"
    H_OT = np.load(path, allow_pickle=True).item()

    mdata["rna"].uns["H_OT"] = H_OT["H_rna"]
    mdata["adt"].uns["H_OT"] = H_OT["H_adt"]

    # Compute the UMAP emebedding
    sc.pp.neighbors(mdata, use_rep="X_mowgli", n_neighbors=30)
    sc.tl.umap(mdata)

    # Compute the silhouette scores
    mdata.obs["silhouette"] = silhouette_samples(
        mdata.obsm["X_mowgli"],
        mdata.obs["rna:celltype"],
    )

    # Plot the cell-wise silhouette score on the UMAP
    sc.pl.umap(mdata, color=["silhouette", "rna:celltype"], show=False)
    plt.savefig(cfg.figure_path + "opcite_silhouette.pdf", bbox_inches="tight")

    # Focus on subsets of CD8 cells
    mdata = mdata[mdata.obs["rna:celltype"].str.contains("CD8+ T", regex=False)]

    # Plot the marker expressions for each celltype
    fig, ax = plt.subplots(figsize=(10, 8), tight_layout=True)
    sc.pl.stacked_violin(
        mdata["adt"],
        [
            "adt:CD8",
            "adt:KLRG1",
            "adt:CD49f",
            "adt:CD57",
            "adt:CD69",
            "adt:TIGIT",
            "adt:CD45RA",
            "adt:CD45RO",
        ],
        groupby="celltype",
        ax=ax,
        use_raw=False,
        show=False,
    )
    plt.savefig(cfg.figure_path + "opcite_markers.pdf", bbox_inches="tight")

    # Look at TIGIT+ cells across batches
    adata = ad.read_h5ad(
        cfg.data_path
        + "OPCITE/GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad"
    )

    # Recover batches that contain TIGIT+CD45RO+ and TIGIT+CD45RA+ cells
    idx = adata.obs["cell_type"].str.contains("TIGIT+ CD45RA+", regex=False)
    contains_cd45ra = adata[idx].obs["Samplename"].unique().tolist()
    idx = adata.obs["cell_type"].str.contains("TIGIT+ CD45RO+", regex=False)
    contains_cd45ro = adata[idx].obs["Samplename"].unique().tolist()
    batches = np.intersect1d(contains_cd45ra, contains_cd45ro)

    for batch in batches:
        # Keep only site one donor 1 becuase of batch effects.
        adata_sub = adata[adata.obs["Samplename"] == batch]

        # Split the variables into RNA and ATAC.
        rna = adata_sub[:, adata_sub.var["feature_types"] == "GEX"].copy()
        adt = adata_sub[:, adata_sub.var["feature_types"] == "ADT"].copy()

        # Rename annotation to follow our convention.
        rna.obs["celltype"] = rna.obs["cell_type"]
        adt.obs["celltype"] = adt.obs["cell_type"]

        # Combine the anndatas into a mudata.
        mdata = mu.MuData({"rna": rna, "adt": adt})

        # Get the TIGIT+ cells.
        tigit_idx = mdata["rna"].obs["celltype"].str.contains("TIGIT+", regex=False)

        fig, axes = plt.subplots(1, 3, figsize=(15, 5), tight_layout=True)

        # CLR preprocessing
        mu.prot.pp.clr(mdata["adt"])

        sc.pl.scatter(
            mdata["adt"][tigit_idx],
            x="CD45RO",
            y="CD45RA",
            color="cell_type",
            size=100,
            alpha=0.7,
            use_raw=False,
            show=False,
            ax=axes[0],
            legend_loc="lower left",
            title=f"Processed counts in {batch.replace('_cite', '').replace('_donor', ' donor ').replace('site', 'site ')}",
        )

        sc.pl.stacked_violin(
            mdata["adt"][tigit_idx],
            ["CD8", "KLRG1", "TIGIT", "CD45RA", "CD45RO"],
            groupby="celltype",
            ax=axes[2],
            use_raw=False,
            show=False,
        )

        mdata["adt"].X = mdata["adt"].layers["counts"]

        sc.pl.scatter(
            mdata["adt"][tigit_idx],
            x="CD45RO",
            y="CD45RA",
            color="cell_type",
            size=100,
            alpha=0.7,
            use_raw=False,
            ax=axes[1],
            show=False,
            legend_loc="lower left",
            title=f"Raw counts in {batch.replace('_cite', '').replace('_donor', ' donor ').replace('site', 'site ')}",
        )
        axes[1].set_yscale("log")
        axes[1].set_xscale("log")

        plt.savefig(cfg.figure_path + f"opcite_{batch}.pdf", bbox_inches="tight")


if __name__ == "__main__":
    my_app()
