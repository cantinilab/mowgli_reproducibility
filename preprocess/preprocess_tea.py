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
    import scanpy as sc
    from rich.console import Console

    ################################## PARAMETERS ########################################
    # Number of highly variable features
    n_highly_variable_genes = 2_500
    n_highly_variable_peaks = 15_000

    ###################################### LOAD DATA #####################################

    console = Console()
    with console.status("[bold green]Loading data..."):

        # Load protein data
        adt = sc.read_csv(
            cfg.data_path
            + "TEA/GSM4949911_X061-AP0C1W1_leukopak_perm-cells"
            + "_tea_fulldepth_adt_counts.csv.gz"
        )

        # Remove the "total" feature.
        adt = adt[:, 1:]

        # Load the 10X multiome data: RNA + ATAC.
        mdata = mu.read_10x_h5(
            cfg.data_path
            + "TEA/GSM4949911_X061-AP0C1W1_leukopak_perm-cells_"
            + "tea_fulldepth_cellranger-arc_filtered_feature_bc_matrix.h5"
        )

        # Make index of proteins compatible with 10X multiome.
        adt.obs.index += "-1"

        # Fuse everything.
        mdata = mu.MuData({"rna": mdata["rna"], "atac": mdata["atac"], "adt": adt})

        # Make variable names unique.
        mdata.var_names_make_unique()

        # Keep only cell with all modalities.
        # (proteins in particular has a lot of background).
        mu.pp.intersect_obs(mdata)

        mdata["rna"].layers["counts"] = mdata["rna"].X.copy()
        mdata["atac"].layers["counts"] = mdata["atac"].X.copy()
        mdata["adt"].layers["counts"] = mdata["adt"].X.copy()

        console.log("Data loaded.")

    ################################## RNA PREPROCESSING #################################

    with console.status("[bold green]Preprocessing RNA..."):

        # Compute quality control metrics.
        mdata["rna"].var["mt"] = mdata["rna"].var_names.str.startswith("rna:MT-")
        sc.pp.calculate_qc_metrics(
            mdata["rna"],
            qc_vars=["mt"],
            percent_top=None,
            log1p=False,
            inplace=True,
        )

        # Filter cells based on QC metrics.
        mu.pp.filter_obs(
            mdata["rna"],
            "n_genes_by_counts",
            lambda x: (x >= 500) & (x < 4_500),
        )
        mu.pp.filter_obs(mdata["rna"], "total_counts", lambda x: x < 12_000)
        mu.pp.filter_obs(mdata["rna"], "pct_counts_mt", lambda x: x < 30)

        # Filter genes by keeping only those that are expressed in at least 10 cells.
        mu.pp.filter_var(mdata["rna"], "n_cells_by_counts", lambda x: x >= 10)

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

    #################################### ATAC PREPROCESSING ##############################

    with console.status("[bold green]Preprocessing ATAC..."):

        # Split interval column into separate variables.
        chr_range = mdata["atac"].var["interval"].str.split(":", expand=True)
        mdata["atac"].var[["Chromosome", "Range"]] = chr_range
        start_end = mdata["atac"].var["Range"].str.split("-", expand=True).astype(int)
        mdata["atac"].var[["Start", "End"]] = start_end
        mdata["atac"].var["Length"] = (
            mdata["atac"].var["End"] - mdata["atac"].var["Start"]
        )

        # Compute QC metrics.
        sc.pp.calculate_qc_metrics(
            mdata["atac"],
            percent_top=None,
            log1p=False,
            inplace=True,
        )

        # Filter peaks based on number of cells where they are present.
        mu.pp.filter_var(mdata["atac"], "n_cells_by_counts", lambda x: x < 4_000)
        mu.pp.filter_var(mdata["atac"], "mean_counts", lambda x: x < 1)
        mu.pp.filter_var(mdata["atac"], "total_counts", lambda x: x < 1e4)
        mu.pp.filter_var(mdata["atac"], "Length", lambda x: x < 5_000)

        # Filter cells based on QC metrics.
        mu.pp.filter_obs(
            mdata["atac"],
            "n_genes_by_counts",
            lambda x: (x >= 1_000) & (x <= 15_000),
        )
        mu.pp.filter_obs(
            mdata["atac"],
            "total_counts",
            lambda x: (x >= 1_000) & (x <= 50_000),
        )

        # Perform per-cell normalization.
        sc.pp.normalize_total(mdata["atac"], target_sum=1e4)
        sc.pp.log1p(mdata["atac"])

        # Only keep highly variable peaks
        sc.pp.highly_variable_genes(
            mdata["atac"],
            n_top_genes=n_highly_variable_peaks,
            subset=True,
            flavor="seurat",
        )

    #################################### ADT PREPROCESSING ###############################

    with console.status("[bold green]Preprocessing ADT..."):

        # Perform per-cell normalization.
        mu.prot.pp.clr(mdata["adt"])
        mdata["adt"].var["highly_variable"] = True

    ################################ SAVE PREPROCESSED DATA ##############################

    with console.status("[bold green]Saving preprocessed data..."):

        # Make sure that signal is not represented as a sparse matrix.
        try:
            mdata["rna"].X = np.array(mdata["rna"].X.todense())
        except:
            pass

        try:
            mdata["atac"].X = np.array(mdata["atac"].X.todense())
        except:
            pass

        try:
            mdata["adt"].X = np.array(mdata["adt"].X.todense())
        except:
            pass

        mu.pp.intersect_obs(mdata)

        # Empty annotation
        mdata["rna"].obs["celltype"] = "none"
        mdata["atac"].obs["celltype"] = "none"
        mdata["adt"].obs["celltype"] = "none"

        # Make variable names unique.
        mdata.var_names_make_unique()

        # Remove this because it prevents from saving.
        del mdata["atac"].uns["atac"]["peak_annotation"]

        # Write the preprocessed data.
        mu.write_h5mu(
            cfg.data_path + "TEA/tea_preprocessed.h5mu.gz",
            mdata,
            compression="gzip",
        )

        console.log("Preprocessed data saved!")


if __name__ == "__main__":
    my_app()
