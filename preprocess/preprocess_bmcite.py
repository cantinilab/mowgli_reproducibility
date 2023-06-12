import hydra
from omegaconf import DictConfig


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def my_app(cfg: DictConfig) -> None:

    ###################################### IMPORTS #######################################

    import warnings

    warnings.simplefilter(action="ignore", category=FutureWarning)
    warnings.simplefilter(action="ignore", category=RuntimeWarning)

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

        # Load the original RNA data.
        rna = sc.read(cfg.data_path + "BMCITE/bm_rna.h5ad")

        # Load the original ADT data.
        adt = sc.read(cfg.data_path + "BMCITE/bm_adt.h5ad")

        # Rename annotation
        rna.obs["celltype"] = rna.obs["celltype.l2"]
        adt.obs["celltype"] = rna.obs["celltype.l2"]

        # Make a MuData object.
        mdata = mu.MuData({"rna": rna, "adt": adt})

        console.log("Metadata loaded.")

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
        mu.pp.filter_obs(
            mdata["rna"],
            "n_genes_by_counts",
            lambda x: (x >= 200) & (x < 2_500),
        )
        mu.pp.filter_obs(
            mdata["rna"],
            "total_counts",
            lambda x: (x >= 800) & (x < 3_500),
        )
        mu.pp.filter_obs(
            mdata["rna"],
            "pct_counts_mt",
            lambda x: x < 5,
        )

        # Filter genes by keeping only those that are expressed in at least 10 cells.
        mu.pp.filter_var(
            mdata["rna"],
            "n_cells_by_counts",
            lambda x: x >= 10,
        )

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

        # Make sure all observations are paired.
        mu.pp.intersect_obs(mdata)

        # Write the preprocessed data.
        mu.write_h5mu(
            cfg.data_path + "BMCITE/bmcite_preprocessed.h5mu.gz",
            mdata,
            compression="gzip",
        )

        # Write the preprocessed RNA data.
        mdata["rna"].write_h5ad(cfg.data_path + "BMCITE/bmcite_preprocessed_rna.h5ad")

        # Write the preprocessed ADT data.
        mdata["adt"].write_h5ad(cfg.data_path + "BMCITE/bmcite_preprocessed_adt.h5ad")

        console.log("Preprocessed data saved!")


if __name__ == "__main__":
    my_app()
