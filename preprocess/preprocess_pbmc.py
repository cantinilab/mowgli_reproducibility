import hydra
from omegaconf import DictConfig


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def my_app(cfg: DictConfig) -> None:

    ######################################## IMPORTS #####################################

    import warnings

    warnings.simplefilter(action="ignore", category=FutureWarning)
    warnings.simplefilter(action="ignore", category=RuntimeWarning)

    # Load libraries.
    import muon as mu
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from rich.console import Console

    ################################## PARAMETERS ########################################
    # Number of highly variable features
    n_highly_variable_genes = 2_500
    n_highly_variable_peaks = 5_000

    ###################################### LOAD DATA #####################################

    console = Console()
    with console.status("[bold green]Loading data..."):

        # Load the original data.
        mdata = mu.read_10x_mtx(
            cfg.data_path + "10X_PBMC_10k/filtered_feature_bc_matrix",
            extended=False,
        )

        console.log("Data loaded.")

        # Add the cell type annotation.
        metadata = pd.read_csv(cfg.data_path + "10X_PBMC_10k/sample_metadata.csv")
        for i in range(metadata.shape[0]):
            if metadata["barcode"][i] in mdata.obs.index:
                mdata.obs.loc[metadata["barcode"][i], "celltype"] = metadata[
                    "celltype"
                ][i]

        # Remove cells without annotation.
        mu.pp.filter_obs(mdata, "celltype", lambda x: pd.notna(x))

        # Copy the cell type annotation to the rna and atac modalities.
        mdata["rna"].obs["celltype"] = mdata.obs["celltype"]
        mdata["atac"].obs["celltype"] = mdata.obs["celltype"]

        mdata["rna"].layers["counts"] = mdata["rna"].X.copy()
        mdata["atac"].layers["counts"] = mdata["atac"].X.copy()

        console.log("Metadata loaded.")

    ################################## RNA PREPROCESSING #################################

    with console.status("[bold green]Preprocessing RNA..."):

        mdata["rna"].layers["counts"] = mdata["rna"].X.copy()
        mdata["atac"].layers["counts"] = mdata["atac"].X.copy()

        # Compute quality control metrics.
        mdata["rna"].var["mt"] = mdata["rna"].var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(
            mdata["rna"], qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )

        # Filter cells based on QC metrics.
        mu.pp.filter_obs(
            mdata["rna"], "n_genes_by_counts", lambda x: (x >= 500) & (x < 4_500)
        )
        mu.pp.filter_obs(mdata["rna"], "total_counts", lambda x: x < 12_000)
        mu.pp.filter_obs(mdata["rna"], "pct_counts_mt", lambda x: x < 20)

        # Filter genes by keeping only those that are expressed in at least 10 cells.
        mu.pp.filter_var(mdata["rna"], "n_cells_by_counts", lambda x: x >= 10)

        # Perform per-cell normalization.
        sc.pp.normalize_total(mdata["rna"], target_sum=1e4)

        # Log-transform the counts.
        sc.pp.log1p(mdata["rna"])

        # Only keep the highly variable genes
        sc.pp.highly_variable_genes(
            mdata["rna"],
            n_top_genes=n_highly_variable_genes,
            subset=True,
            flavor="seurat",
        )

        console.log("RNA preprocessed.")

    ################################## ATAC PREPROCESSING ################################

    with console.status("[bold green]Preprocessing ATAC..."):

        # Compute QC metrics.
        sc.pp.calculate_qc_metrics(
            mdata["atac"], percent_top=None, log1p=False, inplace=True
        )

        # Filter peaks based on number of cells where they are present.
        mu.pp.filter_var(mdata["atac"], "n_cells_by_counts", lambda x: x >= 10)

        # Filter cells based on QC metrics.
        mu.pp.filter_obs(
            mdata["atac"],
            "n_genes_by_counts",
            lambda x: (x >= 2_000) & (x <= 15_000),
        )
        mu.pp.filter_obs(
            mdata["atac"],
            "total_counts",
            lambda x: (x >= 4_000) & (x <= 40_000),
        )

        # Perform per-cell normalization.
        sc.pp.normalize_total(mdata["atac"], target_sum=1e4)

        # Log-transform the counts.
        sc.pp.log1p(mdata["atac"])

        # Only keep highly variable peaks
        sc.pp.highly_variable_genes(
            mdata["atac"],
            n_top_genes=n_highly_variable_peaks,
            subset=True,
            flavor="seurat",
        )

        console.log("ATAC preprocessed.")

    ################################ SAVE PREPROCESSED DATA ##############################

    with console.status("[bold green]Saving data..."):

        # Make sure that signal is not represented as a sparse matrix.
        try:
            mdata["rna"].X = np.array(mdata["rna"].X.todense())
        except:
            pass

        try:
            mdata["atac"].X = np.array(mdata["atac"].X.todense())
        except:
            pass

        # Make sure that all observations are paired.
        mu.pp.intersect_obs(mdata)

        # Write the preprocessed data.
        mu.write_h5mu(
            cfg.data_path + "10X_PBMC_10k/pbmc_preprocessed_with_raw.h5mu.gz",
            mdata,
            compression="gzip",
        )

        console.log("Preprocessed data saved!")


if __name__ == "__main__":
    my_app()
