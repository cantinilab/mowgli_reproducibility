import hydra
from omegaconf import DictConfig


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def my_app(cfg: DictConfig) -> None:

    # Imports.
    import matplotlib.pyplot as plt
    import muon as mu
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from tueplots import axes as tue_axes
    from tueplots import cycler as tue_cycler
    from tueplots import fonts as tue_fonts
    from tueplots.constants.color import palettes as tue_palettes

    import mofax

    sc.set_figure_params(vector_friendly=True, dpi_save=300)

    # Configure the plots.
    plt.rcParams.update({"figure.dpi": 80})
    plt.rcParams.update(tue_axes.spines(left=True, right=False, top=False, bottom=True))
    plt.rcParams.update(tue_axes.grid())
    plt.rcParams.update(tue_cycler.cycler(color=tue_palettes.high_contrast))
    plt.rcParams.update(tue_axes.legend(shadow=False, frameon=False, fancybox=False))
    plt.rcParams.update(tue_fonts.neurips2021_tex(family="sans-serif"))

    # Define the data paths for different datasets.
    data_path = {
        "Liu": cfg.data_path + "Liu/liu_preprocessed.h5mu.gz",
        "sim1": cfg.data_path + "Liu/liu_simulated_1.h5mu.gz",
        "sim2": cfg.data_path + "Liu/liu_simulated_2.h5mu.gz",
        "sim3": cfg.data_path + "Liu/liu_simulated_3.h5mu.gz",
        "sim4": cfg.data_path + "Liu/liu_simulated_4.h5mu.gz",
        "sim5": cfg.data_path + "Liu/liu_simulated_5.h5mu.gz",
        "sim6": cfg.data_path + "Liu/liu_simulated_6.h5mu.gz",
        "10X PBMC": cfg.data_path + "10X_PBMC_10k/pbmc_preprocessed.h5mu.gz",
        "Open Problems Multiome": cfg.data_path
        + "OP_multiome/opmultiome_preprocessed.h5mu.gz",
        "Open Problems CITE-seq": cfg.data_path + "OPCITE/opcite_preprocessed.h5mu.gz",
        "Bone Marrow CITE-seq": cfg.data_path + "BMCITE/bmcite_preprocessed.h5mu.gz",
    }

    # Define the paths to Mowgli's embeddings.
    liu_suffix = "_mowgli_cosine_5_0_1_rna_0_01_atac_0_1_adt_0_01_0_001.npy"
    r_suffix = "_mowgli_cosine_50_0_05_rna_0_01_atac_0_1_adt_0_01_0_001.npy"
    mowgli_path = {
        "Liu": cfg.data_path + "Liu/liu" + liu_suffix,
        "sim1": cfg.data_path + "Liu/liu_simulated_1" + liu_suffix,
        "sim2": cfg.data_path + "Liu/liu_simulated_2" + liu_suffix,
        "sim3": cfg.data_path + "Liu/liu_simulated_3" + liu_suffix,
        "sim4": cfg.data_path + "Liu/liu_simulated_4" + liu_suffix,
        "sim5": cfg.data_path + "Liu/liu_simulated_5" + liu_suffix,
        "sim6": cfg.data_path + "Liu/liu_simulated_6" + liu_suffix,
        "10X PBMC": cfg.data_path + "10X_PBMC_10k/pbmc_preprocessed" + r_suffix,
        "Open Problems Multiome": cfg.data_path + "OP_multiome/opmultiome" + r_suffix,
        "Open Problems CITE-seq": cfg.data_path + "OPCITE/opcite" + r_suffix,
        "Bone Marrow CITE-seq": cfg.data_path + "BMCITE/bmcite" + r_suffix,
    }

    # Define the paths to MOFA+'s embeddings.
    mofa_path = {
        "Liu": cfg.data_path + "Liu/liu_mofa_5.hdf5",
        "sim1": cfg.data_path + "Liu/liu_simulated_1_mofa_5.hdf5",
        "sim2": cfg.data_path + "Liu/liu_simulated_2_mofa_5.hdf5",
        "sim3": cfg.data_path + "Liu/liu_simulated_3_mofa_5.hdf5",
        "sim4": cfg.data_path + "Liu/liu_simulated_4_mofa_5.hdf5",
        "sim5": cfg.data_path + "Liu/liu_simulated_5_mofa_5.hdf5",
        "sim6": cfg.data_path + "Liu/liu_simulated_6_mofa_5.hdf5",
        "10X PBMC": cfg.data_path + "10X_PBMC_10k/pbmc_preprocessed_mofa_15.hdf5",
        "Open Problems Multiome": cfg.data_path + "OP_multiome/opmultiome_mofa_15.hdf5",
        "Open Problems CITE-seq": cfg.data_path + "OPCITE/opcite_mofa_15.hdf5",
        "Bone Marrow CITE-seq": cfg.data_path + "BMCITE/bmcite_mofa_15.hdf5",
    }

    # Define the paths to NMF's embeddings.
    nmf_path = {
        "Liu": cfg.data_path + "Liu/liu_nmf_5.npy",
        "sim1": cfg.data_path + "Liu/liu_simulated_1_nmf_5.npy",
        "sim2": cfg.data_path + "Liu/liu_simulated_2_nmf_5.npy",
        "sim3": cfg.data_path + "Liu/liu_simulated_3_nmf_5.npy",
        "sim4": cfg.data_path + "Liu/liu_simulated_4_nmf_5.npy",
        "sim5": cfg.data_path + "Liu/liu_simulated_5_nmf_5.npy",
        "sim6": cfg.data_path + "Liu/liu_simulated_6_nmf_5.npy",
        "10X PBMC": cfg.data_path + "10X_PBMC_10k/pbmc_preprocessed_nmf_30.npy",
        "Open Problems Multiome": cfg.data_path + "OP_multiome/opmultiome_nmf_30.npy",
        "Open Problems CITE-seq": cfg.data_path + "OPCITE/opcite_nmf_30.npy",
        "Bone Marrow CITE-seq": cfg.data_path + "BMCITE/bmcite_nmf_30.npy",
    }

    # Define the paths to Multigrate's embeddings.
    multigrate_path = {
        "Liu": cfg.data_path + "Liu/liu_multigrate_15.npy",
        "sim1": cfg.data_path + "Liu/liu_simulated_1_multigrate_15.npy",
        "sim2": cfg.data_path + "Liu/liu_simulated_2_multigrate_15.npy",
        "sim3": cfg.data_path + "Liu/liu_simulated_3_multigrate_15.npy",
        "sim4": cfg.data_path + "Liu/liu_simulated_4_multigrate_15.npy",
        "sim5": cfg.data_path + "Liu/liu_simulated_5_multigrate_15.npy",
        "sim6": cfg.data_path + "Liu/liu_simulated_6_multigrate_15.npy",
        "10X PBMC": cfg.data_path + "10X_PBMC_10k/pbmc_preprocessed_multigrate_50.npy",
        "Open Problems Multiome": cfg.data_path
        + "OP_multiome/opmultiome_multigrate_50.npy",
        "Open Problems CITE-seq": cfg.data_path + "OPCITE/opcite_multigrate_50.npy",
        "Bone Marrow CITE-seq": cfg.data_path + "BMCITE/bmcite_multigrate_50.npy",
    }

    # Define the paths to Multigrate's embeddings.
    cobolt_path = {
        "Liu": cfg.data_path + "Liu/liu_cobolt_5.npy",
        "sim1": cfg.data_path + "Liu/liu_simulated_1_cobolt_50.npy",
        "sim2": cfg.data_path + "Liu/liu_simulated_2_cobolt_50.npy",
        "sim3": cfg.data_path + "Liu/liu_simulated_3_cobolt_5.npy",
        "sim4": cfg.data_path + "Liu/liu_simulated_4_cobolt_5.npy",
        "sim5": cfg.data_path + "Liu/liu_simulated_5_cobolt_5.npy",
        "sim6": cfg.data_path + "Liu/liu_simulated_6_cobolt_5.npy",
        "10X PBMC": cfg.data_path + "10X_PBMC_10k/pbmc_preprocessed_cobolt_50.npy",
        "Open Problems Multiome": cfg.data_path
        + "OP_multiome/opmultiome_cobolt_50.npy",
        "Open Problems CITE-seq": cfg.data_path + "OPCITE/opcite_cobolt_50.npy",
        "Bone Marrow CITE-seq": cfg.data_path + "BMCITE/bmcite_cobolt_50.npy",
    }

    # Define the paths to Seurat's UMAP embeddings.
    seurat_path = {
        "Liu": cfg.data_path + "Liu/liu_seurat_umap.csv",
        "sim1": cfg.data_path + "Liu/liu_simulated_1_seurat_umap.csv",
        "sim2": cfg.data_path + "Liu/liu_simulated_2_seurat_umap.csv",
        "sim3": cfg.data_path + "Liu/liu_simulated_3_seurat_umap.csv",
        "sim4": cfg.data_path + "Liu/liu_simulated_4_seurat_umap.csv",
        "sim5": cfg.data_path + "Liu/liu_simulated_5_seurat_umap.csv",
        "sim6": cfg.data_path + "Liu/liu_simulated_6_seurat_umap.csv",
        "10X PBMC": cfg.data_path + "10X_PBMC_10k/pbmc_seurat_umap.csv",
        "Open Problems Multiome": cfg.data_path
        + "OP_multiome/opmultiome_seurat_umap.csv",
        "Open Problems CITE-seq": cfg.data_path + "OPCITE/opcite_seurat_umap.csv",
        "Bone Marrow CITE-seq": cfg.data_path + "BMCITE/bmcite_seurat_umap.csv",
    }

    # Compute and plot the UMAP embeddings
    def plot_umaps(datasets, path, s=15):

        # Define the subplots.
        fig = plt.figure(constrained_layout=True, figsize=(13, 10))
        axes = fig.subplots(len(datasets), 6)

        # Iterate over the datasets.
        for i, dataset in enumerate(datasets):
            print("Dataset: ", dataset)

            # Load the dataset.
            mdata = mu.read_h5mu(data_path[dataset])

            # Iterate over the methods.
            for j in range(6):

                # Plot Mowgli in the 1st column.
                if j == 0:
                    print("> Mowgli")
                    X_mowgli = np.load(mowgli_path[dataset], allow_pickle=True).item()[
                        "W"
                    ]
                    mdata.obsm["X_mowgli"] = X_mowgli
                    mdata.uns = {}
                    sc.pp.neighbors(mdata, use_rep="X_mowgli", n_neighbors=20)
                    sc.tl.umap(mdata)

                # Plot MOFA+ in the 2nd column.
                elif j == 1:
                    print("> MOFA+")
                    mofa_object = mofax.mofa_model(mofa_path[dataset])
                    mdata.obsm["X_mofa"] = mofa_object.get_factors()
                    mdata.uns = {}
                    sc.pp.neighbors(mdata, use_rep="X_mofa", n_neighbors=20)
                    sc.tl.umap(mdata)

                # Plot NMF in the 3nd column.
                elif j == 2:
                    print("> NMF")
                    mdata.obsm["X_nmf"] = np.load(nmf_path[dataset])
                    mdata.uns = {}
                    sc.pp.neighbors(mdata, use_rep="X_nmf", n_neighbors=20)
                    sc.tl.umap(mdata)

                # Plot Seurat in the 4th column.
                elif j == 3:
                    print("> Seurat")
                    mdata.obsm["X_umap"] = pd.read_csv(
                        seurat_path[dataset], index_col=0
                    ).to_numpy()
                    mdata.uns = {}

                # Plot Multigrate in the 5th column.
                elif j == 4:
                    print("> Multigrate")
                    mdata.obsm["X_multigrate"] = np.load(multigrate_path[dataset])
                    mdata.uns = {}
                    sc.pp.neighbors(mdata, use_rep="X_multigrate", n_neighbors=20)
                    sc.tl.umap(mdata)

                # Plot Cobolt in the Tth column.
                elif j == 5:
                    print("> Cobolt")
                    mdata.obsm["X_cobolt"] = np.load(cobolt_path[dataset])
                    mdata.uns = {}
                    sc.pp.neighbors(mdata, use_rep="X_cobolt", n_neighbors=20)
                    sc.tl.umap(mdata)

                # Make the plot.
                sc.pl.umap(
                    mdata,
                    ax=axes[i, j],
                    color="rna:celltype",
                    s=s,
                    alpha=0.7,
                    title="",
                    show=False,
                )

                # Decorate the plots.
                axes[i, j].set_facecolor((0, 0, 0, 0))
                axes[i, j].set(xlabel=None)
                axes[i, j].set(ylabel=None)
                axes[i, j].get_legend().remove()

        # Add the method names as column titles.
        axes[0, 0].set_title("Mowgli")
        axes[0, 1].set_title("MOFA+")
        axes[0, 2].set_title("NMF")
        axes[0, 3].set_title("Seurat")
        axes[0, 4].set_title("Multigrate")
        axes[0, 5].set_title("Cobolt")

        plt.savefig(path)

    # Make the plot for simulated data.
    plot_umaps(
        ["sim1", "sim2", "sim3", "sim4", "sim5", "sim6"],
        cfg.figure_path + "simulated_umaps.pdf",
        s=50,
    )

    # Make the plot for real data.
    plot_umaps(
        [
            "Liu",
            "10X PBMC",
            "Open Problems Multiome",
            "Open Problems CITE-seq",
            "Bone Marrow CITE-seq",
        ],
        cfg.figure_path + "real_umaps.pdf",
    )


if __name__ == "__main__":
    my_app()
