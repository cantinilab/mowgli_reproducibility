import hydra
from omegaconf import DictConfig


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def my_app(cfg: DictConfig) -> None:

    ######################################## IMPORTS #####################################

    import warnings

    warnings.simplefilter(action="ignore", category=FutureWarning)
    warnings.simplefilter(action="ignore", category=RuntimeWarning)

    import muon as mu
    import numpy as np
    import scanpy as sc
    from rich.console import Console

    #################################### PARAMETERS ######################################
    # Number of highly variable features
    n_highly_variable_genes = 1500
    n_highly_variable_peaks = 1500

    ##################################### LOAD DATA ######################################

    console = Console()
    with console.status("[bold green]Loading data..."):

        # Load the original data.
        mdata = mu.read_h5mu(cfg.data_path + "Liu/liu_original.h5mu.gz")

        mdata["rna"].layers["counts"] = mdata["rna"].X.copy()
        mdata["atac"].layers["counts"] = mdata["atac"].X.copy()

        console.log("Raw loaded.")

    ################################## RNA PREPROCESSING #################################

    with console.status("[bold green]Preprocessing RNA..."):

        # Query Ensembl IDs for mitochondrial genes.
        res = sc.queries.mitochondrial_genes("hsapiens", attrname="ensembl_gene_id")
        mito_ensembl_ids = res["ensembl_gene_id"].tolist()

        # Compute quality control metrics.
        mdata["rna"].var["mt"] = [x in mito_ensembl_ids for x in mdata["rna"].var_names]
        sc.pp.calculate_qc_metrics(
            mdata["rna"],
            qc_vars=["mt"],
            percent_top=None,
            log1p=False,
            inplace=True,
        )

        # Filter genes by keeping only those that are expressed in at least 15 cells.
        mu.pp.filter_var(
            mdata["rna"],
            "n_cells_by_counts",
            lambda x: x >= 15,
        )

        # Perform per-cell normalization.
        sc.pp.normalize_total(mdata["rna"])

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
            mdata["atac"],
            percent_top=None,
            log1p=False,
            inplace=True,
        )

        # Filter peaks based on number of cells where they are present.
        mu.pp.filter_var(mdata["atac"], "n_cells_by_counts", lambda x: x >= 15)

        # Filter peaks based on total counts.
        mu.pp.filter_var(mdata["atac"], "total_counts", lambda x: x >= 15)

        # Perform per-cell normalization.
        sc.pp.normalize_total(mdata["atac"])

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

        mu.pp.intersect_obs(mdata)

        # Write the preprocessed data.
        mu.write_h5mu(
            cfg.data_path + "Liu/liu_preprocessed.h5mu.gz", mdata, compression="gzip"
        )

        console.log("Preprocessed data saved!")

        rna_sparsity = np.sum(mdata["rna"].X == 0) / np.size(mdata["rna"].X)
        console.log(f"RNA sparsity: {rna_sparsity}")

        atac_sparsity = np.sum(mdata["atac"].X == 0) / np.size(mdata["atac"].X)
        console.log(f"ATAC sparsity: {atac_sparsity}")

    ################# SIMULATED DATA 1 : mix HCT and HeLa in RNA ##################

    with console.status("[bold green] Simulated 1..."):

        # Load the preprocessed data.
        mdata = mu.read_h5mu(cfg.data_path + "Liu/liu_preprocessed.h5mu.gz")

        # Tamper with the RNA signal
        idx_to_replace = mdata["rna"].obs["celltype"] == "HCT"
        choices = np.where(mdata["rna"].obs["celltype"] == "Hela")[0]
        size = idx_to_replace.sum()
        idx = np.random.choice(choices, size=size)
        mdata["rna"].X[idx_to_replace] = mdata["rna"].X[idx]

        # Filter variables so that none of them are empty
        mdata["atac"].X = np.array(mdata["atac"].X)
        sc.pp.filter_genes(mdata["rna"], min_cells=1)
        sc.pp.filter_genes(mdata["atac"], min_cells=1)

        # Write simulated data
        mdata.write_h5mu(
            cfg.data_path + "Liu/liu_simulated_1.h5mu.gz", compression="gzip"
        )

        console.log("Simulated data saved!")

    ################### SIMULATED DATA 2 : Mix in both omics ######################

    with console.status("[bold green] Simulated 2..."):

        # Load the preprocessed data.
        mdata = mu.read_h5mu(cfg.data_path + "Liu/liu_preprocessed.h5mu.gz")

        # Tamper with the RNA signal.
        idx_to_replace = mdata["rna"].obs["celltype"] == "HCT"
        choices = np.where(mdata["rna"].obs["celltype"] == "Hela")[0]
        size = idx_to_replace.sum()
        idx = np.random.choice(choices, size=size)
        mdata["rna"].X[idx_to_replace] = mdata["rna"].X[idx]

        # Tamper with the ATAC signal.
        idx_to_replace = mdata["atac"].obs["celltype"] == "K562"
        choices = np.where(mdata["atac"].obs["celltype"] == "Hela")[0]
        size = idx_to_replace.sum()
        idx = np.random.choice(choices, size=size)
        mdata["atac"].X[idx_to_replace] = mdata["atac"].X[idx]

        # Filter variables so that none of them are empty
        mdata["atac"].X = np.array(mdata["atac"].X)
        sc.pp.filter_genes(mdata["rna"], min_cells=1)
        sc.pp.filter_genes(mdata["atac"], min_cells=1)

        # Write simulated data
        mdata.write_h5mu(
            cfg.data_path + "Liu/liu_simulated_2.h5mu.gz", compression="gzip"
        )

        console.log("Simulated data saved!")

    #################### SIMULATED DATA 3 : Rare population #######################

    with console.status("[bold green] Simulated 3..."):

        # Load the preprocessed data.
        mdata = mu.read_h5mu(cfg.data_path + "Liu/liu_preprocessed.h5mu.gz")

        # Select the cells to keep.
        idx_hela = np.where(mdata.obs["rna:celltype"] == "Hela")[0]
        idx_hct = np.where(mdata.obs["rna:celltype"] == "HCT")[0]
        idx_k562 = np.where(mdata.obs["rna:celltype"] == "K562")[0]
        random_hela = np.random.choice(idx_hela, size=10, replace=False)
        idx = np.concatenate((random_hela, idx_hct, idx_k562))

        # And only keep these cells.
        mu.pp.filter_obs(mdata, mdata.obs.index[idx])

        # Filter variables so that none of them are empty
        mdata["atac"].X = np.array(mdata["atac"].X)
        sc.pp.filter_genes(mdata["rna"], min_cells=1)
        sc.pp.filter_genes(mdata["atac"], min_cells=1)

        # Write simulated data
        mdata.write_h5mu(
            cfg.data_path + "Liu/liu_simulated_3.h5mu.gz", compression="gzip"
        )

        console.log("Simulated data saved!")

    ################### SIMULATED DATA 4 : Dropout noise ######################

    with console.status("[bold green] Simulated 4..."):

        # Load the preprocessed data.
        mdata = mu.read_h5mu(cfg.data_path + "Liu/liu_preprocessed.h5mu.gz")

        # Add dropout noise to the RNA signal.
        n = mdata["rna"].n_obs
        m = mdata["rna"].var.highly_variable.sum()
        p_dropout = 0.5
        mask = np.random.choice([0, 1], size=(n, m), p=[p_dropout, 1 - p_dropout])
        mdata["rna"].X[:, mdata["rna"].var.highly_variable] *= mask

        rna_sparsity = np.sum(mdata["rna"].X == 0) / np.size(mdata["rna"].X)
        console.log(f"rna sparsity: {rna_sparsity}")

        # Add dropout noise to the ATAC signal.
        n = mdata["atac"].n_obs
        m = mdata["atac"].var.highly_variable.sum()
        p_dropout = 0.5
        mask = np.random.choice([0, 1], size=(n, m), p=[p_dropout, 1 - p_dropout])
        mdata["atac"].X[:, mdata["atac"].var.highly_variable] *= mask

        atac_sparsity = np.sum(mdata["atac"].X == 0) / np.size(mdata["atac"].X)
        console.log(f"atac sparsity: {atac_sparsity}")

        # Filter variables so that none of them are empty
        mdata["atac"].X = np.array(mdata["atac"].X)
        sc.pp.filter_genes(mdata["rna"], min_cells=1)
        sc.pp.filter_genes(mdata["atac"], min_cells=1)

        # Write simulated data
        mdata.write_h5mu(
            cfg.data_path + "Liu/liu_simulated_4.h5mu.gz", compression="gzip"
        )

        console.log("Simulated data saved!")

    ################### SIMULATED DATA 5 : Dropout noise ######################

    with console.status("[bold green] Simulated 5..."):

        # Load the preprocessed data.
        mdata = mu.read_h5mu(cfg.data_path + "Liu/liu_preprocessed.h5mu.gz")

        # Add dropout noise to the RNA signal.
        n = mdata["rna"].n_obs
        m = mdata["rna"].var.highly_variable.sum()
        p_dropout = 0.7
        mask = np.random.choice([0, 1], size=(n, m), p=[p_dropout, 1 - p_dropout])
        mdata["rna"].X[:, mdata["rna"].var.highly_variable] *= mask

        rna_sparsity = np.sum(mdata["rna"].X == 0) / np.size(mdata["rna"].X)
        console.log(f"rna sparsity: {rna_sparsity}")

        # Add dropout noise to the ATAC signal.
        n = mdata["atac"].n_obs
        m = mdata["atac"].var.highly_variable.sum()
        p_dropout = 0.7
        mask = np.random.choice([0, 1], size=(n, m), p=[p_dropout, 1 - p_dropout])
        mdata["atac"].X[:, mdata["atac"].var.highly_variable] *= mask

        atac_sparsity = np.sum(mdata["atac"].X == 0) / np.size(mdata["atac"].X)
        console.log(f"atac sparsity: {atac_sparsity}")

        # Filter variables so that none of them are empty
        mdata["atac"].X = np.array(mdata["atac"].X)
        sc.pp.filter_genes(mdata["rna"], min_cells=1)
        sc.pp.filter_genes(mdata["atac"], min_cells=1)

        # Write simulated data
        mdata.write_h5mu(
            cfg.data_path + "Liu/liu_simulated_5.h5mu.gz", compression="gzip"
        )

        console.log("Simulated data saved!")

    ################### SIMULATED DATA 6 : Dropout noise ######################

    with console.status("[bold green] Simulated 6..."):

        # Load the preprocessed data.
        mdata = mu.read_h5mu(cfg.data_path + "Liu/liu_preprocessed.h5mu.gz")

        # Add dropout noise to the RNA signal.
        n = mdata["rna"].n_obs
        m = mdata["rna"].var.highly_variable.sum()
        p_dropout = 0.9
        mask = np.random.choice([0, 1], size=(n, m), p=[p_dropout, 1 - p_dropout])
        mdata["rna"].X[:, mdata["rna"].var.highly_variable] *= mask

        rna_sparsity = np.sum(mdata["rna"].X == 0) / np.size(mdata["rna"].X)
        console.log(f"rna sparsity: {rna_sparsity}")

        # Add dropout noise to the ATAC signal.
        n = mdata["atac"].n_obs
        m = mdata["atac"].var.highly_variable.sum()
        p_dropout = 0.9
        mask = np.random.choice([0, 1], size=(n, m), p=[p_dropout, 1 - p_dropout])
        mdata["atac"].X[:, mdata["atac"].var.highly_variable] *= mask

        atac_sparsity = np.sum(mdata["atac"].X == 0) / np.size(mdata["atac"].X)
        console.log(f"atac sparsity: {atac_sparsity}")

        # Filter variables so that none of them are empty
        mdata["atac"].X = np.array(mdata["atac"].X)
        sc.pp.filter_genes(mdata["rna"], min_cells=1)
        sc.pp.filter_genes(mdata["atac"], min_cells=1)

        # Write simulated data
        mdata.write_h5mu(
            cfg.data_path + "Liu/liu_simulated_6.h5mu.gz", compression="gzip"
        )

        console.log("Simulated data saved!")


if __name__ == "__main__":
    my_app()
