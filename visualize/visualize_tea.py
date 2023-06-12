import hydra
from omegaconf import DictConfig


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def my_app(cfg: DictConfig) -> None:

    import anndata as ad
    import matplotlib.pyplot as plt
    import muon as mu
    import numpy as np
    import pandas as pd
    import plotly.graph_objects as go
    import scanpy as sc
    import seaborn as sns
    from matplotlib.ticker import FormatStrFormatter

    import mofax

    sc.set_figure_params(vector_friendly=True, dpi_save=300)

    ##########################################################################################
    #################################### Loading the data ####################################
    ##########################################################################################
    # Load the raw RNA signal.
    rna_all_genes = mu.read_10x_h5(
        cfg.data_path
        + "TEA/GSM4949911_X061-AP0C1W1_leukopak_perm-cells_"
        + "tea_fulldepth_cellranger-arc_filtered_feature_bc_matrix.h5",
        extended=False,
    )["rna"]

    # Compute quality control metrics and perform basic preprocessing.
    rna_all_genes.var["mt"] = rna_all_genes.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        rna_all_genes, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    mu.pp.filter_obs(
        rna_all_genes, "n_genes_by_counts", lambda x: (x >= 500) & (x < 4_500)
    )
    mu.pp.filter_obs(rna_all_genes, "total_counts", lambda x: x < 12_000)
    mu.pp.filter_obs(rna_all_genes, "pct_counts_mt", lambda x: x < 30)
    mu.pp.filter_var(rna_all_genes, "n_cells_by_counts", lambda x: x >= 10)
    sc.pp.normalize_total(
        rna_all_genes, target_sum=1e4
    )  # Perform per-cell normalization.
    sc.pp.log1p(rna_all_genes)  # Log-transform the counts.

    # Load the preprocessed data.
    mdata = mu.read_h5mu(cfg.data_path + "TEA/tea_preprocessed.h5mu.gz")

    # This is needed somehow.
    mdata.uns = {}

    # Subset rna_all_genes to the cells in mdata.
    rna_all_genes = rna_all_genes[mdata["rna"].obs_names]
    rna_all_genes.var_names_make_unique()

    ##########################################################################################
    ######################################## Load MOFA #######################################
    ##########################################################################################

    # Load the MOFA model.
    mofa_model = mofax.mofa_model(cfg.data_path + "TEA/tea_mofa_15.hdf5")
    mdata.obsm["X_mofa_15"] = mofa_model.get_factors()

    # Load MOFA+'s loadings.
    H_mofa = {
        "H_rna": mofa_model.get_weights("rna"),
        "H_atac": mofa_model.get_weights("atac"),
        "H_adt": mofa_model.get_weights("adt"),
    }

    mofa_model = mofax.mofa_model(cfg.data_path + "TEA/tea_mofa_30.hdf5")
    mdata.obsm["X_mofa_30"] = mofa_model.get_factors()

    mofa_model = mofax.mofa_model(cfg.data_path + "TEA/tea_mofa_50.hdf5")
    mdata.obsm["X_mofa_50"] = mofa_model.get_factors()

    ##########################################################################################
    ###################################### Load Mowgli #######################################
    ##########################################################################################

    # Load the Mowgli model.
    mowgli_path = (
        cfg.data_path
        + "TEA/tea_mowgli_cosine_15_0_05_rna_0_01_atac_0_1_adt_0_01_0_001.npy"
    )
    mdata.obsm["X_mowgli_15"] = np.load(mowgli_path, allow_pickle=True).item()["W"]

    mowgli_path = (
        cfg.data_path
        + "TEA/tea_mowgli_cosine_30_0_05_rna_0_01_atac_0_1_adt_0_01_0_001.npy"
    )
    mdata.obsm["X_mowgli_30"] = np.load(mowgli_path, allow_pickle=True).item()["W"]

    mowgli_path = (
        cfg.data_path
        + "TEA/tea_mowgli_cosine_50_0_05_rna_0_01_atac_0_1_adt_0_01_0_01.npy"
    )
    mdata.obsm["X_mowgli_w_0_01"] = np.load(mowgli_path, allow_pickle=True).item()["W"]

    mowgli_path = (
        cfg.data_path
        + "TEA/tea_mowgli_cosine_50_0_05_rna_0_01_atac_0_1_adt_0_01_0_0001.npy"
    )
    mdata.obsm["X_mowgli_w_0_0001"] = np.load(mowgli_path, allow_pickle=True).item()[
        "W"
    ]

    mowgli_path = (
        cfg.w_path + "tea_mowgli_cosine_50_0_05_rna_0_01_atac_0_1_adt_0_01_0_001.npy"
    )
    mdata.obsm["X_mowgli_50"] = np.load(mowgli_path, allow_pickle=True).item()["W"]

    # Load Mowgli's loadings.
    H_mowgli = np.load(
        cfg.h_path + "tea_mowgli_cosine_50_0_05_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
        allow_pickle=True,
    ).item()

    ##########################################################################################
    ######################################## Load NMF ########################################
    ##########################################################################################

    # Load the NMF model.
    nmf = np.load(cfg.data_path + "TEA/tea_nmf_15.npy", allow_pickle=True).item()
    mdata.obsm["X_nmf_15"] = nmf["W"]

    nmf = np.load(cfg.data_path + "TEA/tea_nmf_30.npy", allow_pickle=True).item()
    mdata.obsm["X_nmf_30"] = nmf["W"]

    # Load NMF's loadings.
    H_nmf = {
        "H_rna": nmf["H_rna"].T,
        "H_atac": nmf["H_atac"].T,
        "H_adt": nmf["H_adt"].T,
    }

    nmf = np.load(cfg.data_path + "TEA/tea_nmf_50.npy", allow_pickle=True).item()
    mdata.obsm["X_nmf_50"] = nmf["W"]

    ##########################################################################################
    ##################################### Define markers #####################################
    ##########################################################################################

    # Define protein markers for cell types.
    adt_markers = {
        "Eryth": ["CD71"],
        "B": ["CD19", "CD21", "IgD", "IgM"],
        "Mono": ["CD14", "CD141", "CD11b", "CD172a", "CD11c"],
        "NK": ["CD16", "CD56"],
        "T": ["CD3"],
        "MAIT": ["TCR-Va7.2"],
        "CD8 T": ["CD8a"],
        "CD4 T": ["CD4"],
    }

    # Define a flat list of adt markers.
    adt_markers_flat = [m for markers in adt_markers.values() for m in markers]

    # Define RNA markers for cell types.
    rna_markers = {
        "Eryth": ["HBD", "HBM", "TRIM58"],
        "B": ["RALGPS2", "MS4A1", "BANK1", "IGHM"],
        "Mono": ["CTSS", "FCN1", "LYZ", "PSAP", "S100A9"],
        "NK": ["KLRD1", "KLRF1", "CST7", "GZMB", "NKG7", "GNLY"],
        "T": ["CD3D", "CD3G", "TRAC"],
        "MAIT": ["KLRB1", "GZMK", "SLC4A10", "GZMA"],
        "CD8 T": ["CD8A", "CD8B", "LINC02446"],
        "CD4 T": ["CD4"],
    }

    # Make sure that RNA markers are present in our data.
    varnames = rna_all_genes.var_names
    for celltype in rna_markers:
        rna_markers[celltype] = [g for g in rna_markers[celltype] if g in varnames]
    rna_markers_flat = [m for markers in rna_markers.values() for m in markers]

    # Define Mowgli marker factors for cell types.
    mowgli_factor_markers = {
        "Eryth": ["7"],
        "B": ["33"],
        "Mono": ["32"],
        "NK": ["2"],
        "MAIT": ["9"],
        "CD8": ["16", "49"],
        "CD4": ["8"],
    }

    # Define MOFA+ positive marker factors for cell types.
    mofa_pos_factor_markers = {
        "B": ["0"],
        "Mono": ["1"],
        "NK": ["3", "2"],
        "MAIT": ["4"],
        "CD8": ["6", "9"],
    }

    # Define MOFA+ negative marker factors for cell types.
    mofa_neg_factor_markers = {
        "Eryth": ["5", "10"],
        "Mono": ["2"],
        "MAIT": ["9"],
        "CD4": ["0"],
    }

    nmf_factor_markers = {
        "Eryth": ["9"],
        "B": ["28"],
        "Mono": ["12"],
        "NK": ["6", "16"],
        "MAIT": ["0"],
        "CD8": ["8", "21"],
        "CD4": ["2", "22"],
    }

    # Load the azimuth annotation.
    azimuth_pred = pd.read_csv(
        cfg.data_path + "TEA/azimuth_pred.tsv", sep="\t", index_col=0
    )
    mdata.obs[azimuth_pred.columns] = azimuth_pred.loc[mdata.obs_names]

    ##########################################################################################
    ######################################### Annotate #######################################
    ##########################################################################################

    def annotate(
        method="mowgli",
        cluster_order=None,
        method_cluster_names=None,
        resolution=0.2,
        cmap="Purples",
    ):

        print("Annotating with", method)

        # Compute neighbors UMAP embedding and Leiden clustering for the method.
        nb_key = f"{method}_neighbors"
        cl_key = f"leiden_{method}"
        min_dist = {
            "mowgli_50": 0.8,
            "mofa_15": 0.5,
            "nmf_30": 0.6,
        }
        sc.pp.neighbors(mdata, n_neighbors=20, key_added=nb_key, use_rep=f"X_{method}")
        sc.tl.umap(mdata, neighbors_key=nb_key, min_dist=min_dist[method])
        sc.tl.leiden(
            mdata, resolution=resolution, key_added=cl_key, neighbors_key=nb_key
        )

        # Keyword params.
        dotplot_kwds = {
            "groupby": cl_key,
            "mean_only_expressed": True,
            "size_title": "Fraction of cells\nin cluster (%)",
            "categories_order": cluster_order,
            "show": False,
            "cmap": cmap,
        }
        umap_kwds = {
            "alpha": 0.7,
            "legend_fontoutline": 2,
            "frameon": False,
            "show": False,
        }

        # Make the UMAP plot, colored by cluster.
        sc.pl.umap(
            mdata,
            color=cl_key,
            legend_loc="on data",
            title=f"Leiden clustering of {method} embedding",
            **umap_kwds,
        )
        plt.savefig(cfg.figure_path + f"{method}_tea_umap.pdf", bbox_inches="tight")
        plt.close()

        print("Computed UMAP and Leiden clustering for", method)
        print("There are", len(mdata.obs[cl_key].unique()), "clusters")

        # Make a dotplot with ADT values for each cluster.
        mdata["adt"].obs[cl_key] = mdata.obs[cl_key]
        mdata["adt"].var_names = mdata["adt"].var_names.str.replace("adt:", "")
        axes = sc.pl.dotplot(
            mdata["adt"],
            adt_markers,
            title=f"{method}: counts for marker proteins in each cluster",
            colorbar_title="Mean ADT count",
            expression_cutoff=0.5,
            **dotplot_kwds,
        )
        axes["mainplot_ax"].set_ylabel("Cluster #")
        axes["mainplot_ax"].set_xlabel("Marker proteins")
        plt.savefig(
            cfg.figure_path + f"{method}_tea_leiden_adt.pdf", bbox_inches="tight"
        )
        plt.close()

        # Make a dotplot with RNA values for each cluster.
        rna_all_genes.obs[cl_key] = mdata.obs[cl_key]
        axes = sc.pl.dotplot(
            rna_all_genes,
            rna_markers,
            title=f"{method}: counts for marker genes in each cluster",
            colorbar_title="Mean gene counts",
            **dotplot_kwds,
        )
        axes["mainplot_ax"].set_ylabel("Cluster #")
        axes["mainplot_ax"].set_xlabel("Marker genes")
        plt.savefig(
            cfg.figure_path + f"{method}_tea_leiden_rna.pdf", bbox_inches="tight"
        )
        plt.close()

        # Annotate the method's embedding.
        mdata.obs[f"annotation_{method}"] = mdata.obs["predicted.celltype.l2"]
        if method_cluster_names:
            codes = mdata.obs[cl_key].cat.codes
            mdata.obs[f"annotation_{method}"] = [method_cluster_names[c] for c in codes]

        # Make a UMAP plot of the method's embedding, colored by annotation.
        sc.pl.umap(
            mdata,
            color=f"annotation_{method}",
            title=f"Annotated {method} embedding",
            legend_fontweight="normal",
            **umap_kwds,
        )
        plt.savefig(
            cfg.figure_path + f"{method}_tea_umap_annotated.pdf", bbox_inches="tight"
        )
        plt.close()

    annotate(
        method="mowgli_50",
        cluster_order=["8", "0", "5", "6", "7", "3", "1", "2", "4"],
        method_cluster_names={
            0: "B cells",
            1: "CD4 T cells",
            2: "CD4 T cells",
            3: "CD8 T cells",
            4: "CD4 T cells",
            5: "Monocytes",
            6: "NK cells",
            7: "MAIT T cells",
            8: "Erythroid cells",
        },
    )

    annotate(
        method="mofa_15",
        cluster_order=["9", "1", "8", "3", "5", "7", "6", "4", "0", "2"],
        method_cluster_names={
            0: "CD4 T cells",
            1: "B cells",
            2: "CD4 T cells",
            3: "Monocytes",
            4: "CD8 T cells",
            5: "NK cells",
            6: "CD8 T cells",
            7: "MAIT T cells",
            8: "B cells",
            9: "Erythroid cells",
        },
    )

    annotate(
        method="nmf_30",
        cluster_order=[
            "11",
            "2",
            "9",
            "10",
            "3",
            "6",
            "8",
            "5",
            "7",
            "0",
            "1",
            "4",
        ],
        method_cluster_names={
            0: "CD4 T cells",
            1: "CD4 T cells",
            2: "B cells",
            3: "Monocytes",
            4: "CD4 T cells",
            5: "CD8 T cells",
            6: "NK cells",
            7: "CD8 T cells",
            8: "MAIT T cells",
            9: "B cells",
            10: "B cells",
            11: "Erythroid cells",
        },
        resolution=1.0,
    )

    ####################### Sankey diagram of the cell type annotations ######################

    annot_azimuth = mdata.obs["predicted.celltype.l1"]
    annot_mofa = mdata.obs["annotation_mofa_15"]
    annot_mowgli = mdata.obs["annotation_mowgli_50"]
    annot_nmf = mdata.obs["annotation_nmf_30"]

    agreement = annot_mofa == annot_mowgli
    agreement &= annot_nmf == annot_mowgli
    print("Agreement between Mowgli, NMF and MOFA+:", np.mean(agreement))

    categories_mowgli = np.sort(annot_mowgli.cat.categories).tolist()
    categories_azimuth = np.sort(annot_azimuth.cat.categories).tolist()
    categories_mofa = np.sort(annot_mofa.cat.categories).tolist()
    categories_nmf = np.sort(annot_nmf.cat.categories).tolist()

    categories = (
        categories_mowgli + categories_azimuth + categories_mofa + categories_nmf
    )
    source, target, value = [], [], []

    for i, a in enumerate(categories_mowgli):
        for j, b in enumerate(categories_azimuth):
            source.append(i)
            target.append(len(categories_mowgli) + j)
            value.append(np.sum((annot_mowgli == a) & (annot_azimuth == b)))

    for i, a in enumerate(categories_azimuth):
        for j, b in enumerate(categories_mofa):
            source.append(len(categories_mowgli) + i)
            target.append(len(categories_mowgli) + len(categories_azimuth) + j)
            value.append(np.sum((annot_azimuth == a) & (annot_mofa == b)))

    for i, a in enumerate(categories_mofa):
        for j, b in enumerate(categories_nmf):
            source.append(len(categories_mowgli) + len(categories_azimuth) + i)
            target.append(
                len(categories_mowgli)
                + len(categories_azimuth)
                + len(categories_mofa)
                + j
            )
            value.append(np.sum((annot_mofa == a) & (annot_nmf == b)))

    # Colors for Mowgli
    colors = [
        "#002626",
        "#0E4749",
        "#E6E075",
        "#F07C42",
        "#C2A370",
        "#142A58",
        "#C87E7E",
    ]

    # Colors for Azimuth
    colors += [
        "#002626",
        "#0E4749",
        "#E6E075",
        "#F36C6C",
        "#142A58",
        "#C87E7E",
        "#F07C42",
        "#C2A370",
    ]

    # Colors for MOFA+
    colors += [
        "#002626",
        "#0E4749",
        "#E6E075",
        "#F07C42",
        "#C2A370",
        "#142A58",
        "#C87E7E",
    ]

    # Colors for NMF
    colors += [
        "#002626",
        "#0E4749",
        "#E6E075",
        "#F07C42",
        "#C2A370",
        "#142A58",
        "#C87E7E",
    ]

    line = {"color": "black", "width": 0.5}
    node = {
        "pad": 10,
        "thickness": 20,
        "line": line,
        "label": categories,
        "color": colors,
    }
    link = {"source": source, "target": target, "value": value}

    fig = go.Figure(data=[go.Sankey(node=node, link=link)])
    fig.update_layout(font_size=10, width=500, height=500)
    fig.write_image(cfg.figure_path + "sankey_mowgli_azimuth_mofa_nmf.pdf")

    ##########################################################################################
    ############################ Interpret the methods' dimensions ###########################
    ##########################################################################################

    def interpret(
        method,
        method_factor_markers,
        adt_dict,
        rna_dict,
        categories_order=[
            "Erythroid cells",
            "B cells",
            "Monocytes",
            "NK cells",
            "MAIT T cells",
            "CD8 T cells",
            "CD4 T cells",
        ],
        hint="",
        cmap="Purples",
    ):
        method_factor_markers_flat = [
            m for l in method_factor_markers.values() for m in l
        ]

        # Make a dotplot of weights for the method's factors across clusters.
        method_embedding = ad.AnnData(mdata.obsm[f"X_{method}"])
        method_embedding.obs_names = mdata.obs_names
        method_embedding.obs[f"annotation_{method}"] = mdata.obs[f"annotation_{method}"]
        axes = sc.pl.dotplot(
            method_embedding,
            method_factor_markers,
            groupby=f"annotation_{method}",
            categories_order=categories_order,
            expression_cutoff=1e-2,
            mean_only_expressed=True,
            title=f"Weights of {method}'s factors",
            colorbar_title="Mean weight",
            size_title="Fraction of cells\nin cluster (%)",
            cmap="Purples",
            show=False,
        )
        axes["mainplot_ax"].set_ylabel("Cluster")
        axes["mainplot_ax"].set_xlabel("Factor #")
        axes["mainplot_ax"].set_xticklabels(
            axes["mainplot_ax"].get_xticklabels(), rotation=0
        )
        plt.savefig(
            cfg.figure_path + f"{method}{hint}_tea_factors.pdf", bbox_inches="tight"
        )
        plt.close()

        # Make a matrixplot of ADT weights accross the method's factors.
        adata = ad.AnnData(adt_dict)
        adata.X /= adata.X.std(axis=1, keepdims=True)
        adata.X[adata.X < 0] = 0
        adata.obs_names = mdata["adt"].var_names
        adata.obs["adt"] = pd.Categorical(adata.obs_names)
        adata = adata[adt_markers_flat, method_factor_markers_flat]
        sc.pl.matrixplot(
            adata,
            method_factor_markers,
            groupby="adt",
            cmap=cmap,
            categories_order=adt_markers_flat,
            title="Proteins weights in the method's top factors",
            colorbar_title="ADT weight",
            show=False,
        )
        plt.savefig(
            cfg.figure_path + f"{method}{hint}_tea_factors_adt.pdf", bbox_inches="tight"
        )
        plt.close()

        adata = ad.AnnData(rna_dict)
        adata.X /= adata.X.std(axis=1, keepdims=True)
        adata.X[adata.X < 0] = 0
        adata.obs_names = mdata["rna"].var_names.str.replace("rna:", "")
        adata.obs["rna"] = pd.Categorical(adata.obs_names)
        genes = [g for g in rna_markers_flat if g in adata.obs_names]
        adata = adata[genes, method_factor_markers_flat]
        sc.pl.matrixplot(
            adata,
            method_factor_markers,
            groupby="rna",
            cmap=cmap,
            categories_order=genes,
            title="Gene weights in the method's top factors",
            colorbar_title="Gene weight",
            show=False,
        )
        plt.savefig(
            cfg.figure_path + f"{method}{hint}_tea_factors_rna.pdf", bbox_inches="tight"
        )
        plt.close()

    interpret(
        method="mowgli_50",
        method_factor_markers=mowgli_factor_markers,
        adt_dict=H_mowgli["H_adt"].copy(),
        rna_dict=H_mowgli["H_rna"].copy(),
        cmap="Purples",
    )

    interpret(
        method="mofa_15",
        hint="_pos",
        method_factor_markers=mofa_pos_factor_markers,
        adt_dict=H_mofa["H_adt"].copy(),
        rna_dict=H_mofa["H_rna"].copy(),
        cmap="Blues",
    )

    interpret(
        method="mofa_15",
        hint="_neg",
        method_factor_markers=mofa_neg_factor_markers,
        adt_dict=-H_mofa["H_adt"].copy(),
        rna_dict=-H_mofa["H_rna"].copy(),
        cmap="Reds",
    )

    interpret(
        method="nmf_30",
        method_factor_markers=nmf_factor_markers,
        adt_dict=H_nmf["H_adt"].copy(),
        rna_dict=H_nmf["H_rna"].copy(),
        cmap="Oranges",
    )

    ##########################################################################################
    ###################################### Bubble plots ######################################
    ##########################################################################################

    def plot_bubbles(
        clusters,
        embedding_key,
        annotation_key,
        colors,
        save_to="tea_factors_bubble.pdf",
    ):

        fig, axes = plt.subplots(
            len(embedding_key),
            len(clusters),
            figsize=(3.5 * len(clusters), 3 * len(embedding_key)),
            constrained_layout=True,
        )

        for row_id, method in enumerate(embedding_key.keys()):

            for i, cluster in enumerate(clusters):

                df = []

                annot_key = annotation_key[method]

                embedding_pos = ad.AnnData(mdata.obsm[embedding_key[method]])
                embedding_pos.obs_names = mdata.obs_names
                embedding_pos.obs[annot_key] = mdata.obs[annot_key]
                embedding_pos.X[embedding_pos.X < 0] = 0

                embedding_neg = ad.AnnData(-mdata.obsm[embedding_key[method]])
                embedding_neg.obs_names = mdata.obs_names
                embedding_neg.obs[annot_key] = mdata.obs[annot_key]
                embedding_neg.X[embedding_neg.X < 0] = 0

                for factor in embedding_pos.var_names:

                    idx_cl = embedding_pos.obs[annot_key] == cluster
                    embed_ravel = embedding_pos[idx_cl, factor].X.ravel()
                    mean_cluster = embed_ravel.mean()
                    mean_other = embedding_pos[~idx_cl, factor].X.ravel().mean()
                    prop_expressed = float(np.mean(embed_ravel > 1e-3))
                    if prop_expressed > 0:
                        df.append(
                            {
                                "factor": factor,
                                "Sign of weights": "Positive",
                                "mean_cluster": float(mean_cluster),
                                "mean_other": float(mean_other),
                                "Proportion": float(prop_expressed),
                            }
                        )

                for factor in embedding_neg.var_names:
                    idx_cl = embedding_neg.obs[annot_key] == cluster
                    embed_ravel = embedding_neg[idx_cl, factor].X.ravel()
                    mean_cluster = embed_ravel.mean()
                    mean_other = embedding_neg[~idx_cl, factor].X.ravel().mean()
                    prop_expressed = float(np.mean(embed_ravel > 1e-3))
                    if prop_expressed > 0:
                        df.append(
                            {
                                "factor": factor,
                                "Sign of weights": "Negative",
                                "mean_cluster": float(mean_cluster),
                                "mean_other": float(mean_other),
                                "Proportion": float(prop_expressed),
                            }
                        )

                df = pd.DataFrame(df)

                max_val = df[["mean_other", "mean_cluster"]].to_numpy().max()
                padding = 0.1 * max_val
                axes[row_id, i].plot(
                    [0, max_val],
                    [0, max_val],
                    linewidth=2,
                    linestyle="--",
                    color="grey",
                )

                sns.scatterplot(
                    data=df,
                    x="mean_other",
                    y="mean_cluster",
                    hue="Sign of weights",
                    palette=colors[method],
                    size="Proportion",
                    ax=axes[row_id, i],
                    sizes=(500 * df["Proportion"].min(), 500 * df["Proportion"].max()),
                    alpha=0.7,
                )

                for j in range(len(df)):
                    score = df["mean_cluster"].iloc[j] - df["mean_other"].iloc[j]
                    score /= df["mean_cluster"].max()
                    if score > 0.25:
                        axes[row_id, i].text(
                            df["mean_other"].iloc[j],
                            df["mean_cluster"].iloc[j],
                            df["factor"].iloc[j],
                            horizontalalignment="center",
                            verticalalignment="center",
                            size="medium",
                            color="white",
                            weight="semibold",
                            alpha=0.8,
                        )

                scores = (df["mean_cluster"] - df["mean_other"]) / df[
                    "mean_cluster"
                ].max()
                for j, score in enumerate(scores):
                    if score > 0.6:
                        axes[row_id, i].text(
                            df["mean_other"].iloc[j] + padding,
                            df["mean_cluster"].iloc[j],
                            f"score = {score:.2f}",
                            horizontalalignment="left",
                            verticalalignment="center",
                            size="medium",
                            color="black",
                            weight="semibold",
                            alpha=0.8,
                        )

                axes[row_id, i].set_ylabel("Mean weight in cell type")
                axes[row_id, i].set_xlabel("Mean weight other cell types")

                if i < len(clusters) - 1:
                    axes[row_id, i].get_legend().remove()
                else:
                    axes[row_id, i].legend(
                        loc="upper left",
                        bbox_to_anchor=(1.05, 1),
                        frameon=False,
                        title="Prop. of cells in cluster with\nabsolute weights > 1e-3",
                    )
                    for lh in axes[row_id, i].get_legend().legendHandles:
                        lh.set_alpha(0.7)

                axes[row_id, i].spines.right.set_visible(False)
                axes[row_id, i].spines.top.set_visible(False)
                axes[row_id, i].yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
                axes[row_id, i].xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
                axes[row_id, i].set_xlim(-padding, max_val + padding)
                axes[row_id, i].set_ylim(-padding, max_val + padding)

        for i, cluster in enumerate(clusters):
            axes[0, i].set_title(cluster)

        plt.savefig(cfg.figure_path + save_to, bbox_inches="tight")
        plt.close()

    clusters = ["Erythroid cells", "B cells", "Monocytes", "NK cells"]
    clusters += ["MAIT T cells", "CD8 T cells", "CD4 T cells"]

    plot_bubbles(
        clusters=clusters,
        embedding_key={
            "Mowgli": "X_mowgli_50",
            "MOFA+": "X_mofa_15",
            "NMF": "X_nmf_30",
        },
        annotation_key={
            "Mowgli": "annotation_mowgli_50",
            "MOFA+": "annotation_mofa_15",
            "NMF": "annotation_nmf_30",
        },
        colors={
            "Mowgli": ["purple"],
            "MOFA+": ["tab:blue", "tab:red"],
            "NMF": ["peru"],
        },
    )

    plot_bubbles(
        clusters=clusters,
        embedding_key={
            "Mowgli_0_01": "X_mowgli_w_0_01",
            "Mowgli_0_001": "X_mowgli_50",
            "Mowgli_0_0001": "X_mowgli_w_0_0001",
        },
        annotation_key={
            "Mowgli_0_01": "annotation_mowgli_50",
            "Mowgli_0_001": "annotation_mowgli_50",
            "Mowgli_0_0001": "annotation_mowgli_50",
        },
        colors={
            "Mowgli_0_01": ["purple"],
            "Mowgli_0_001": ["purple"],
            "Mowgli_0_0001": ["purple"],
        },
        save_to="mowgli_w_tea_factors_bubble.pdf",
    )

    plot_bubbles(
        clusters=clusters,
        embedding_key={
            "MOFA+ (15 factors)": "X_mofa_15",
            "MOFA+ (30 factors)": "X_mofa_30",
            "MOFA+ (50 factors)": "X_mofa_50",
        },
        annotation_key={
            "MOFA+ (15 factors)": "annotation_mofa_15",
            "MOFA+ (30 factors)": "annotation_mofa_15",
            "MOFA+ (50 factors)": "annotation_mofa_15",
        },
        colors={
            "MOFA+ (15 factors)": ["tab:blue", "tab:red"],
            "MOFA+ (30 factors)": ["tab:blue", "tab:red"],
            "MOFA+ (50 factors)": ["tab:blue", "tab:red"],
        },
        save_to="mofa_latent_tea_factors_bubble.pdf",
    )

    plot_bubbles(
        clusters=clusters,
        embedding_key={
            "Mowgli (15 factors)": "X_mowgli_15",
            "Mowgli (30 factors)": "X_mowgli_30",
            "Mowgli (50 factors)": "X_mowgli_50",
        },
        annotation_key={
            "Mowgli (15 factors)": "annotation_mowgli_50",
            "Mowgli (30 factors)": "annotation_mowgli_50",
            "Mowgli (50 factors)": "annotation_mowgli_50",
        },
        colors={
            "Mowgli (15 factors)": ["purple"],
            "Mowgli (30 factors)": ["purple"],
            "Mowgli (50 factors)": ["purple"],
        },
        save_to="mowgli_latent_tea_factors_bubble.pdf",
    )

    plot_bubbles(
        clusters=clusters,
        embedding_key={
            "NMF (15 factors)": "X_nmf_15",
            "NMF (30 factors)": "X_nmf_30",
            "NMF (50 factors)": "X_nmf_50",
        },
        annotation_key={
            "NMF (15 factors)": "annotation_nmf_30",
            "NMF (30 factors)": "annotation_nmf_30",
            "NMF (50 factors)": "annotation_nmf_30",
        },
        colors={
            "NMF (15 factors)": ["peru"],
            "NMF (30 factors)": ["peru"],
            "NMF (50 factors)": ["peru"],
        },
        save_to="nmf_latent_tea_factors_bubble.pdf",
    )

    # Describe the marker-sepcific factors

    for i in range(50):
        mdata["adt"].obs[f"mowgli_{i}"] = mdata.obsm["X_mowgli_50"][:, i]

    H = np.load(
        cfg.h_path + "tea_mowgli_cosine_50_0_05_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
        allow_pickle=True,
    ).item()

    for k in H:
        mdata.uns[k] = H[k]

    mdata["adt"].obs[mdata["adt"].var_names] = mdata["adt"].X

    sc.tl.umap(mdata, neighbors_key="mowgli_50_neighbors", min_dist=0.8)
    mdata["adt"].obsm["X_umap"] = mdata.obsm["X_umap"]

    # CD56 bright
    # NK cell that has the phenotype CD56-bright, CD16-negative, and CD84-positive with the function to secrete interferon-gamma but is not cytotoxic. [database_cross_reference: PMID:19796267][database_cross_reference: GO_REF:0000031][database_cross_reference: GOC:pam][database_cross_reference: GOC:add][database_cross_reference: GOC:tfm][database_cross_reference: PMID:22343568]

    fig, axes = plt.subplots(6, 4, figsize=(30, 30))

    # barplot
    H_adt = mdata.uns["H_adt"]
    H_adt /= H_adt.sum(axis=1, keepdims=True)
    xx = H_adt[:, 13]
    yy = mdata["adt"].var_names
    idx = np.argsort(xx)[::-1]
    sns.barplot(x=xx[idx[:20]], y=yy[idx[:20]], ax=axes[0, 0], palette="Blues_r")

    # scatterplot
    sns.scatterplot(
        x=mdata["adt"].obs["CD56"],
        y=mdata["adt"].obs["CD16"],
        hue=mdata["adt"].obs["mowgli_13"],
        s=15,
        alpha=0.7,
        palette="viridis",
        ax=axes[0, 1],
        rasterized=True,
    )

    mdata["adt"].obs["adt:CD56"] = mdata["adt"].obs["CD56"]
    sc.pl.umap(mdata["adt"], color="adt:CD56", ax=axes[0, 2], show=False, alpha=0.7)
    sc.pl.umap(mdata["adt"], color="mowgli_13", ax=axes[0, 3], show=False, alpha=0.7)

    # nonclassical monocyte
    # A patrolling monocyte that is CD14-low and CD16-positive. [database_cross_reference: PMID:20870168][database_cross_reference: GOC:tfm]

    # barplot
    H_adt = mdata.uns["H_adt"]
    H_adt /= H_adt.sum(axis=1, keepdims=True)
    xx = H_adt[:, 34]
    yy = mdata["adt"].var_names
    idx = np.argsort(xx)[::-1]
    sns.barplot(x=xx[idx[:20]], y=yy[idx[:20]], ax=axes[1, 0], palette="Blues_r")

    # scatterplot
    sns.scatterplot(
        x=mdata["adt"].obs["CD14"],
        y=mdata["adt"].obs["CD16"],
        hue=mdata["adt"].obs["mowgli_34"],
        s=15,
        alpha=0.7,
        palette="viridis",
        ax=axes[1, 1],
        rasterized=True,
    )

    mdata["adt"].obs["adt:CD16"] = mdata["adt"].obs["CD16"]
    sc.pl.umap(mdata["adt"], color="adt:CD16", ax=axes[1, 2], show=False, alpha=0.7)
    sc.pl.umap(mdata["adt"], color="mowgli_34", ax=axes[1, 3], show=False, alpha=0.7)

    # classical monocyte
    # A classical monocyte that is CD14-positive, CD16-negative, CD64-positive, CD163-positive. [database_cross_reference: PMID:19689341][database_cross_reference: PMID:15615263][database_cross_reference: PMID:1706877][database_cross_reference: GOC:add][database_cross_reference: GOC:tfm][database_cross_reference: PMID:22343568]
    # expression of CCR2 in both rodents and humans, negative for the lineage markers CD3, CD19, and CD20,

    # barplot
    H_adt = mdata.uns["H_adt"]
    H_adt /= H_adt.sum(axis=1, keepdims=True)
    xx = H_adt[:, 32]
    yy = mdata["adt"].var_names
    idx = np.argsort(xx)[::-1]
    sns.barplot(x=xx[idx[:20]], y=yy[idx[:20]], ax=axes[2, 0], palette="Blues_r")

    # scatterplot
    sns.scatterplot(
        x=mdata["adt"].obs["CD14"],
        y=mdata["adt"].obs["CD16"],
        hue=mdata["adt"].obs["mowgli_32"],
        s=15,
        alpha=0.7,
        palette="viridis",
        ax=axes[2, 1],
        rasterized=True,
    )

    mdata["adt"].obs["adt:CD14"] = mdata["adt"].obs["CD14"]
    sc.pl.umap(mdata["adt"], color="adt:CD14", ax=axes[2, 2], show=False, alpha=0.7)
    sc.pl.umap(mdata["adt"], color="mowgli_32", ax=axes[2, 3], show=False, alpha=0.7)

    # cDC2: A myeloid dendritic cell found in the blood that is CD1c-positive. [database_cross_reference: PMID:20628149][database_cross_reference: PMID:20204387][database_cross_reference: GOC:tfm][database_cross_reference: GOC:dsd]
    # Normally represent 10-20% of peripheral blood mDCs (human). They are also CD281-positive (TLR1), CD282-positive (TLR2), CD283-positive (TLR3), CD284-positive (TLR4), CD285-positive (TLR5), CD286-positive (TLR6), CD288-positive (TLR8), and CD290-positive (TLR10) [PMID:20204387]. Upon TLR stimulation, these cells were potent producers of CXCL8 (IL-8), while producing little TNF-alpha.

    # barplot
    H_adt = mdata.uns["H_adt"]
    H_adt /= H_adt.sum(axis=1, keepdims=True)
    xx = H_adt[:, 15]
    yy = mdata["adt"].var_names
    idx = np.argsort(xx)[::-1]
    sns.barplot(x=xx[idx[:20]], y=yy[idx[:20]], ax=axes[3, 0], palette="Blues_r")

    # scatterplot
    sns.scatterplot(
        x=mdata["adt"].obs["mowgli_15"],
        y=mdata["adt"].obs["FceRI"],
        ax=axes[3, 1],
        s=10,
        alpha=0.7,
        rasterized=True,
    )

    mdata["adt"].obs["adt:FceRI"] = mdata["adt"].obs["FceRI"]
    sc.pl.umap(mdata["adt"], color="adt:FceRI", ax=axes[3, 2], show=False, alpha=0.7)
    sc.pl.umap(mdata["adt"], color="mowgli_15", ax=axes[3, 3], show=False, alpha=0.7)

    # pDC
    # DC definition: A cell of hematopoietic origin, typically resident in particular tissues, specialized in the uptake, processing, and transport of antigens to lymph nodes for the purpose of stimulating an immune response via T cell activation. These cells are lineage negative (CD3-negative, CD19-negative, CD34-negative, and CD56-negative). [database_cross_reference: GOC:add][database_cross_reference: ISBN:0781735149]
    # A dendritic cell type of distinct morphology, localization, and surface marker expression (CD123-positive) from other dendritic cell types and associated with early stage immune responses, particularly the release of physiologically abundant amounts of type I interferons in response to infection. [database_cross_reference: GOC:add][database_cross_reference: PMCID:PMC538703][database_cross_reference: GOC:dsd][database_cross_reference: PMID:17850486][database_cross_reference: PMID:15549123][database_cross_reference: PMID:20304825][database_cross_reference: PMID:17332250][database_cross_reference: PMCID:PMC2118448]

    # barplot
    H_adt = mdata.uns["H_adt"]
    H_adt /= H_adt.sum(axis=1, keepdims=True)
    xx = H_adt[:, 3]
    yy = mdata["adt"].var_names
    idx = np.argsort(xx)[::-1]
    sns.barplot(x=xx[idx[:20]], y=yy[idx[:20]], ax=axes[4, 0], palette="Blues_r")

    # scatterplot
    sns.scatterplot(
        x=mdata["adt"].obs["mowgli_3"],
        y=mdata["adt"].obs["CD304"],
        ax=axes[4, 1],
        s=10,
        alpha=0.7,
        rasterized=True,
    )
    mdata["adt"].obs["adt:CD304"] = mdata["adt"].obs["CD304"]
    sc.pl.umap(mdata["adt"], color="adt:CD304", ax=axes[4, 2], show=False, alpha=0.7)
    sc.pl.umap(mdata["adt"], color="mowgli_3", ax=axes[4, 3], show=False, alpha=0.7)

    # barplot
    H_adt = mdata.uns["H_adt"]
    H_adt /= H_adt.sum(axis=1, keepdims=True)
    xx = H_adt[:, 8]
    yy = mdata["adt"].var_names
    idx = np.argsort(xx)[::-1]
    sns.barplot(x=xx[idx[:20]], y=yy[idx[:20]], ax=axes[5, 0], palette="Blues_r")

    # scatterplot
    sns.scatterplot(
        x=mdata["adt"].obs["CD4"],
        y=np.log1p(mdata["adt"].obs["CD45RO"] / mdata["adt"].obs["CD45RA"]),
        hue=mdata["adt"].obs["mowgli_8"],
        s=10,
        alpha=0.9,
        palette="viridis",
        ax=axes[5, 1],
        rasterized=True,
    )
    axes[5, 1].set_ylabel("log(1 + CD45RO / CD45RA)")

    mdata["adt"].obs["adt:CD45RA"] = mdata["adt"].obs["CD45RA"]
    sc.pl.umap(mdata["adt"], color="adt:CD45RA", ax=axes[5, 2], show=False, alpha=0.7)
    sc.pl.umap(mdata["adt"], color="mowgli_8", ax=axes[5, 3], show=False, alpha=0.7)

    plt.savefig(cfg.figure_path + "specific_tea_umap.pdf", bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    my_app()
