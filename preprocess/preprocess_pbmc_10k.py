######################################## IMPORTS #########################################

# Load libraries.
import numpy as np
import pandas as pd
import scanpy as sc
import muon as mu
from rich.console import Console
import os


################################## PARAMETERS ########################################

# Define data path.
data_path = os.path.join(
    "/users/csb/huizing/Documents/PhD/Code/",
    "mowgli_reproducibility/data/10X_PBMC_10k",
)

# Number of highly variable features
n_highly_variable_genes = 2_500
n_highly_variable_peaks = 5_000

console = Console()
with console.status("[bold green]Loading data...") as status:

    ###################################### LOAD DATA #####################################

    # Load the original data.
    mdata = mu.read_10x_mtx(
        os.path.join(data_path, "filtered_feature_bc_matrix"),
        extended=False,
    )

    console.log("Data loaded.")

    # Add the cell type annotation.
    metadata = pd.read_csv(os.path.join(data_path, "sample_metadata.csv"))
    for i in range(metadata.shape[0]):
        if metadata["barcode"][i] in mdata.obs.index:
            mdata.obs.loc[metadata["barcode"][i], "celltype"] = metadata["celltype"][i]

    # Remove cells without annotation.
    mu.pp.filter_obs(mdata, "celltype", lambda x: pd.notna(x))

    # Copy the cell type annotation to the rna and atac modalities.
    mdata["rna"].obs["celltype"] = mdata.obs["celltype"]
    mdata["atac"].obs["celltype"] = mdata.obs["celltype"]

    console.log("Metadata loaded.")

rna, atac = mdata["rna"], mdata["atac"]

with console.status("[bold green]Preprocessing data...") as status:

    ################################## RNA PREPROCESSING #################################

    print(rna.var_names)

    # Compute quality control metrics.
    rna.var["mt"] = rna.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        rna, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    # Filter cells based on QC metrics.
    mu.pp.filter_obs(rna, "n_genes_by_counts", lambda x: (x >= 500) & (x < 4_500))
    mu.pp.filter_obs(rna, "total_counts", lambda x: x < 12_000)
    mu.pp.filter_obs(rna, "pct_counts_mt", lambda x: x < 20)

    # Filter genes by keeping only those that are expressed in at least 10 cells.
    mu.pp.filter_var(rna, "n_cells_by_counts", lambda x: x >= 10)

    # Perform per-cell normalization.
    sc.pp.normalize_total(rna, target_sum=1e4)

    # Log-transform the counts.
    sc.pp.log1p(rna)

    # Only keep the highly variable genes
    sc.pp.highly_variable_genes(
        rna, n_top_genes=n_highly_variable_genes, subset=True, flavor="seurat"
    )

    console.log("RNA preprocessed.")

    ################################## ATAC PREPROCESSING ################################

    # Compute QC metrics.
    sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)

    # Filter peaks based on number of cells where they are present.
    mu.pp.filter_var(atac, "n_cells_by_counts", lambda x: x >= 10)

    # Filter cells based on QC metrics.
    mu.pp.filter_obs(
        atac,
        "n_genes_by_counts",
        lambda x: (x >= 2_000) & (x <= 15_000),
    )
    mu.pp.filter_obs(
        atac,
        "total_counts",
        lambda x: (x >= 4_000) & (x <= 40_000),
    )

    # Perform per-cell normalization.
    sc.pp.normalize_total(atac, target_sum=1e4)

    # Log-transform the counts.
    sc.pp.log1p(atac)

    # Only keep highly variable peaks
    sc.pp.highly_variable_genes(
        atac,
        n_top_genes=n_highly_variable_peaks,
        subset=True,
        flavor="seurat",
    )

    console.log("ATAC preprocessed.")

with console.status("[bold green]Saving data...") as status:

    mdata = mu.MuData({"rna": rna, "atac": atac})

    ################################ SAVE PREPROCESSED DATA ##############################

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
        os.path.join(data_path, "pbmc_preprocessed.h5mu.gz"),
        mdata,
        compression="gzip",
    )
    rna.write_h5ad(
        os.path.join(data_path, "pbmc_preprocessed_rna.h5ad.gz"),
        compression="gzip",
    )
    atac.write_h5ad(
        os.path.join(data_path, "pbmc_preprocessed_atac.h5ad.gz"),
        compression="gzip",
    )

    console.log("Preprocessed data saved!")
