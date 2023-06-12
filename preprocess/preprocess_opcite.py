import hydra
from omegaconf import DictConfig


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def my_app(cfg: DictConfig) -> None:

    #################################### IMPORTS #########################################

    import warnings

    warnings.simplefilter(action="ignore", category=FutureWarning)
    warnings.simplefilter(action="ignore", category=RuntimeWarning)

    import anndata as ad
    import muon as mu
    import numpy as np
    import scanpy as sc
    from rich.console import Console

    ################################## PARAMETERS ########################################
    # Number of highly variable features
    n_highly_variable_genes = 2_500

    ###################################### LOAD DATA #####################################

    console = Console()
    with console.status("[bold green]Loading data..."):

        # Load the original data.
        adata = ad.read_h5ad(
            cfg.data_path
            + "OPCITE/GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad"
        )

        # Keep only site one donor 1 becuase of batch effects.
        adata = adata[adata.obs["Samplename"] == "site1_donor1_cite"]

        # Split the variables into RNA and ATAC.
        rna = adata[:, adata.var["feature_types"] == "GEX"].copy()
        adt = adata[:, adata.var["feature_types"] == "ADT"].copy()

        # Rename annotation to follow our convention.
        rna.obs["celltype"] = rna.obs["cell_type"]
        adt.obs["celltype"] = adt.obs["cell_type"]

        # Combine the anndatas into a mudata.
        mdata = mu.MuData({"rna": rna, "adt": adt})

        console.log("Data loaded.")

    ################################## RNA PREPROCESSING #################################

    with console.status("[bold green]Preprocessing RNA..."):

        # Compute quality control metrics.
        mdata["rna"].var["mt"] = mdata["rna"].var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(
            mdata["rna"],
            qc_vars=["mt"],
            percent_top=None,
            log1p=False,
            inplace=True,
        )

        # Filter cells based on QC metrics.
        mu.pp.filter_obs(mdata["rna"], "n_genes_by_counts", lambda x: x < 5_000)
        mu.pp.filter_obs(mdata["rna"], "total_counts", lambda x: x < 20_000)
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

    ################################## ADT PREPROCESSING #################################

    with console.status("[bold green]Preprocessing ADT..."):

        # Perform Center Log Ratio preprocessing.
        mu.prot.pp.clr(mdata["adt"])

        mdata["adt"].var["highly_variable"] = True

        console.log("ADT preprocessed.")

    ################################ SAVE PREPROCESSED DATA ##############################

    with console.status("[bold green]Saving data..."):

        # Make sure that signal is not represented as a sparse matrix.
        try:
            mdata["rna"].X = np.array(mdata["rna"].X.todense())
        except:
            pass

        try:
            mdata["adt"].X = np.array(mdata["adt"].X.todense())
        except:
            pass

        # Make variable names unique.
        mdata.var_names_make_unique()

        mu.pp.intersect_obs(mdata)

        # Write the preprocessed data.
        mu.write_h5mu(
            cfg.data_path + "OPCITE/opcite_preprocessed.h5mu.gz",
            mdata,
            compression="gzip",
        )

        # Write the preprocessed RNA data.
        mdata["rna"].write_h5ad(cfg.data_path + "OPCITE/opcite_preprocessed_rna.h5ad")

        # Write the preprocessed ADT data.
        mdata["adt"].write_h5ad(cfg.data_path + "OPCITE/opcite_preprocessed_adt.h5ad")

        console.log("Preprocessed data saved!")


if __name__ == "__main__":
    my_app()
