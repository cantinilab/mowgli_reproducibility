##########################################################################################
######################################## Imports #########################################
##########################################################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import muon as mu
import scanpy as sc
import mofax
import anndata as ad
from scipy.stats import pearsonr
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter
from matplotlib.lines import Line2D

##########################################################################################
#################################### Define some paths ###################################
##########################################################################################

# Define the folder containing data.
data_folder = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/"

# Define the folder where to save figures.
fig_folder = "/users/csb/huizing/Documents/PhD/Code/"
fig_folder += "mowgli_reproducibility/visualize/figures/"

# Define the folder for enrichment results.
enrich_folder = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/enrich/"

# Define the folder containing Mowgli's embeddings.
w_folder = "/users/csb/huizing/Documents/PhD/Code/Mowgli/local_analysis/from_jz/w/"

# Define the folder containing Mowgli's loadings.
h_folder = "/users/csb/huizing/Documents/PhD/Code/Mowgli/local_analysis/from_jz/h/"

##########################################################################################
#################################### Loading the data ####################################
##########################################################################################

# Load the raw RNA signal.
rna_all_genes = mu.read_10x_h5(
    data_folder
    + "TEA/GSM4949911_X061-AP0C1W1_leukopak_perm-cells_"
    + "tea_fulldepth_cellranger-arc_filtered_feature_bc_matrix.h5",
    extended=False,
)["rna"]

# Compute quality control metrics and perform basic preprocessing.
rna_all_genes.var["mt"] = rna_all_genes.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    rna_all_genes,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True,
)
mu.pp.filter_obs(rna_all_genes, "n_genes_by_counts", lambda x: (x >= 500) & (x < 4_500))
mu.pp.filter_obs(rna_all_genes, "total_counts", lambda x: x < 12_000)
mu.pp.filter_obs(rna_all_genes, "pct_counts_mt", lambda x: x < 30)
mu.pp.filter_var(rna_all_genes, "n_cells_by_counts", lambda x: x >= 10)
sc.pp.normalize_total(rna_all_genes, target_sum=1e4)  # Perform per-cell normalization.
sc.pp.log1p(rna_all_genes)  # Log-transform the counts.

# Load the preprocessed data.
mdata = mu.read_h5mu(data_folder + "TEA/tea_preprocessed.h5mu.gz")

# Load the MOFA model.
mofa_model = mofax.mofa_model(data_folder + "TEA/tea_mofa_15.hdf5")
mdata.obsm["X_mofa"] = mofa_model.get_factors()

# Load the Mowgli model.
mdata.obsm["X_mowgli"] = np.load(
    w_folder + "tea_mowgli_cosine_50_0_05_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    allow_pickle=True,
).item()["W"]

# This is needed somehow.
mdata.uns = {}

# Subset rna_all_genes to the cells in mdata.
rna_all_genes = rna_all_genes[mdata["rna"].obs_names]
rna_all_genes.var_names_make_unique()

# Load MOFA+'s loadings.
H_mofa = {
    "H_rna": mofa_model.get_weights("rna"),
    "H_atac": mofa_model.get_weights("atac"),
    "H_adt": mofa_model.get_weights("adt"),
}

# Load Mowgli's loadings.
H_mowgli = np.load(
    h_folder + "tea_mowgli_cosine_50_0_05_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    allow_pickle=True,
).item()

##########################################################################################
################################# Load enrichment scores #################################
##########################################################################################

# Define the method used for enrichement.
enrichment_file = "gprofiler"  # or "fgsea" or "enrichr"

# Get fgsea results.
if enrichment_file == "fgsea":

    # Load the fgsea results, rename columns and make dimensions integers.
    enr = pd.read_csv(enrich_folder + "fgsea.csv", index_col=0)
    cols = {"pathway": "native", "col": "dim", "pval": "p_value"}
    enr.rename(columns=cols, inplace=True)
    enr["dim"] = enr["dim"].str.replace("X", "").astype(int)

# Get gprofiler results.
elif enrichment_file == "gprofiler":

    # Load the gprofiler results, rename columns and make dimensions integers.
    enr = pd.read_csv(enrich_folder + "enrichment.csv", index_col=0)
    enr["dim"] = enr["query"].str.split(" ", expand=True)[1].astype(int)

# Get enrichr results.
elif enrichment_file == "enrichr":

    # Load the enrichr results, rename columns and make dimensions integers.
    enr = pd.read_csv(enrich_folder + "enrichr.csv", index_col=0)
    enr.rename(columns={"path_name": "native", "adj_p_val": "p_value"}, inplace=True)

# Add some more columns.
enr["name"] = enr["native"]
enr = enr[enr["p_value"] < 0.05]
enr["minlogp"] = -np.log10(enr["p_value"])

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

# Define a flat list of rna markers.
rna_markers_flat = [m for markers in rna_markers.values() for m in markers]

# Define Mowgli marker factors for cell types.
mowgli_factor_markers = {
    "Eryth": ["7", "30"],
    "B": ["33", "44"],
    "Mono": ["32", "34"],
    "NK": ["2"],
    "MAIT": ["9"],
    "CD8": ["16", "49"],
    "CD4": ["8", "18", "28"],
}

# Define a flat list of Mowgli marker factors.
mowgli_factor_markers_flat = [m for l in mowgli_factor_markers.values() for m in l]

# Add the rest of the factors.
# mowgli_factor_markers[" "] = [
# str(i) for i in range(50) if str(i) not in mowgli_factor_markers_flat
# ]
# mowgli_factor_markers_flat += mowgli_factor_markers[" "]

# Define MOFA+ positive marker factors for cell types.
mofa_pos_factor_markers = {
    "Eryth": ["6"],
    "B": ["0", "11"],
    "Mono": ["1", "8"],
    "NK": ["3"],
    "MAIT": ["4"],
    "CD8": ["9"],
}

# Define a flat list of MOFA+ positive marker factors.
mofa_pos_factor_markers_flat = [m for l in mofa_pos_factor_markers.values() for m in l]

# Add the rest of the factors.
# mofa_pos_factor_markers[" "] = [
# str(i) for i in range(15) if str(i) not in mofa_pos_factor_markers_flat
# ]
# mofa_pos_factor_markers_flat += mofa_pos_factor_markers[" "]

# Define MOFA+ negative marker factors for cell types.
mofa_neg_factor_markers = {
    "Eryth": ["5", "10"],
    "Mono": ["2", "8"],
    "MAIT": ["9", "3"],
}

# Define a flat list of MOFA+ negative marker factors.
mofa_neg_factor_markers_flat = [m for l in mofa_neg_factor_markers.values() for m in l]

# Add the rest of the factors.
# mofa_neg_factor_markers[" "] = [
# str(i) for i in range(15) if str(i) not in mofa_neg_factor_markers_flat
# ]
# mofa_neg_factor_markers_flat += mofa_neg_factor_markers[" "]

# Define gene sets for cell types.
gene_sets = {
    "B gene set": "B Cell CL0000785",
    "Mono gene set": "Monocyte CL0000576",
    "NK gene set": "Natural Killer CL0000623",
    "MAIT gene set": "Mucosal Associated Invariant T CL0000940",
    "CD8 gene set": "CD8 T CL0000625",
    "CD4 gene set": "CD4 T CL0000624",
}

##########################################################################################
##################################### Annotate Mowgli ####################################
##########################################################################################

# Compute neighbors UMAP embedding and Leiden clustering for Mowgli.
sc.pp.neighbors(mdata, n_neighbors=25, key_added="mowgli_neighbors", use_rep="X_mowgli")
sc.tl.umap(mdata, neighbors_key="mowgli_neighbors")
sc.tl.leiden(
    mdata,
    resolution=0.2,
    key_added="leiden_mowgli",
    neighbors_key="mowgli_neighbors",
)

# Make the UMAP plot, colored by cluster.
ax = sc.pl.umap(
    mdata,
    color="leiden_mowgli",
    alpha=0.7,
    legend_loc="on data",
    title="Leiden clustering of Mowgli embedding",
    legend_fontoutline=2,
    frameon=False,
    show=False,
)
plt.savefig(fig_folder + "mowgli_tea_umap.pdf", bbox_inches="tight")
plt.close()

# Make a dotplot with ADT values for each cluster.
mdata["adt"].obs["leiden_mowgli"] = mdata.obs["leiden_mowgli"]
mdata["adt"].var_names = mdata["adt"].var_names.str.replace("adt:", "")
axes = sc.pl.dotplot(
    mdata["adt"],
    adt_markers,
    groupby="leiden_mowgli",
    title="Mowgli: counts for marker proteins in each cluster",
    mean_only_expressed=True,
    colorbar_title="Mean ADT count",
    categories_order=["8", "0", "5", "6", "7", "3", "1", "2", "4"],
    size_title="Fraction of cells\nin cluster (%)",
    expression_cutoff=0.5,
    show=False,
    cmap="Purples",
)
axes["mainplot_ax"].set_ylabel("Cluster #")
axes["mainplot_ax"].set_xlabel("Marker proteins")
plt.savefig(fig_folder + "mowgli_tea_leiden_adt.pdf", bbox_inches="tight")
plt.close()

# Make a dotplot with RNA values for each cluster.
rna_all_genes.obs["leiden_mowgli"] = mdata.obs["leiden_mowgli"]
axes = sc.pl.dotplot(
    rna_all_genes,
    rna_markers,
    groupby="leiden_mowgli",
    title="Mowgli: counts for marker genes in each cluster",
    categories_order=["8", "0", "5", "6", "7", "3", "1", "2", "4"],
    mean_only_expressed=True,
    colorbar_title="Mean gene counts",
    size_title="Fraction of cells\nin cluster (%)",
    show=False,
    cmap="Purples",
)
axes["mainplot_ax"].set_ylabel("Cluster #")
axes["mainplot_ax"].set_xlabel("Marker genes")
plt.savefig(fig_folder + "mowgli_tea_leiden_rna.pdf", bbox_inches="tight")
plt.close()

# Annotate the Mowgli embedding.
mowgli_cluster_names = {
    0: "B cells",
    1: "CD4 T cells",
    2: "CD4 T cells",
    3: "CD8 T cells",
    4: "CD4 T cells",
    5: "Monocytes",
    6: "NK cells",
    7: "MAIT T cells",
    8: "Erythroid cells",
}
codes = mdata.obs["leiden_mowgli"].cat.codes
mdata.obs["annotation_mowgli"] = [mowgli_cluster_names[c] for c in codes]

# Make a UMAP plot of the Mowgli embedding, colored by annotation.
ax = sc.pl.umap(
    mdata,
    color="annotation_mowgli",
    alpha=0.7,
    legend_loc="on data",
    title="Annotated Mowgli embedding",
    legend_fontoutline=2,
    frameon=False,
    legend_fontweight="normal",
    show=False,
)
plt.savefig(fig_folder + "mowgli_tea_umap_annotated.pdf", bbox_inches="tight")
plt.close()

##########################################################################################
############################## Interpret Mowgli's dimensions #############################
##########################################################################################

# Make a dotplot of weights for Mowgli's factors across clusters.
mowgli_embedding = ad.AnnData(mdata.obsm["X_mowgli"])
mowgli_embedding.obs_names = mdata.obs_names
mowgli_embedding.obs["annotation_mowgli"] = mdata.obs["annotation_mowgli"]
axes = sc.pl.dotplot(
    mowgli_embedding,
    mowgli_factor_markers,
    groupby="annotation_mowgli",
    categories_order=[
        "Erythroid cells",
        "B cells",
        "Monocytes",
        "NK cells",
        "MAIT T cells",
        "CD8 T cells",
        "CD4 T cells",
    ],
    expression_cutoff=1e-2,
    mean_only_expressed=True,
    title="Weights of Mowgli's factors",
    colorbar_title="Mean weight",
    size_title="Fraction of cells\nin cluster (%)",
    cmap="Purples",
    show=False,
)
axes["mainplot_ax"].set_ylabel("Cluster")
axes["mainplot_ax"].set_xlabel("Factor #")
axes["mainplot_ax"].set_xticklabels(axes["mainplot_ax"].get_xticklabels(), rotation=0)
plt.savefig(fig_folder + "mowgli_tea_factors.pdf", bbox_inches="tight")
plt.close()

# TODO: Make a bubble plot of Mowgli's weights for each factor, CD4 vs other clusters


# Make a matrixplot of ADT weights accross Mowgli's factors.
adata = ad.AnnData(H_mowgli["H_adt"].copy())
adata.X /= adata.X.std(axis=1, keepdims=True)
adata.obs_names = mdata["adt"].var_names
adata.obs["adt"] = pd.Categorical(adata.obs_names)
adata = adata[adt_markers_flat, mowgli_factor_markers_flat]
sc.pl.matrixplot(
    adata,
    mowgli_factor_markers,
    groupby="adt",
    cmap="Purples",
    categories_order=adt_markers_flat,
    title="Proteins weights in Mowgli's top factors",
    colorbar_title="ADT weight",
    # standard_scale="group",
    show=False,
)
plt.savefig(fig_folder + "mowgli_tea_factors_adt.pdf", bbox_inches="tight")
plt.close()


adata = ad.AnnData(H_mowgli["H_rna"].copy())
adata.X /= adata.X.std(axis=1, keepdims=True)
adata.obs_names = mdata["rna"].var_names.str.replace("rna:", "")
adata.obs["rna"] = pd.Categorical(adata.obs_names)
genes = [g for g in rna_markers_flat if g in adata.obs_names]
adata = adata[genes, mowgli_factor_markers_flat]
sc.pl.matrixplot(
    adata,
    mowgli_factor_markers,
    groupby="rna",
    cmap="Purples",
    categories_order=genes,
    title="Gene weights in Mowgli's top factors",
    colorbar_title="Gene weight",
    # standard_scale="group",
    show=False,
)
plt.savefig(fig_folder + "mowgli_tea_factors_rna.pdf", bbox_inches="tight")
plt.close()

mowgli_pvals = ad.AnnData(np.zeros((len(gene_sets), mowgli_embedding.n_vars)))
mowgli_pvals.obs_names = gene_sets.keys()
mowgli_pvals.var_names = mowgli_embedding.var_names

for i, x in enumerate(gene_sets):
    idx = enr["native"] == gene_sets[x]
    idx &= enr["method"] == "mowgli"
    mowgli_pvals.X[i, enr.loc[idx, "dim"]] = enr.loc[idx, "minlogp"]

mowgli_pvals.obs["gene_set"] = pd.Categorical(mowgli_pvals.obs_names)
sc.pl.matrixplot(
    mowgli_pvals,
    mowgli_factor_markers,
    groupby="gene_set",
    categories_order=gene_sets.keys(),
    cmap="Purples",
    title="Cell type enrichment for the factors",
    colorbar_title=r"$-\log_{10}(p~value)$",
    show=False,
)
plt.savefig(fig_folder + "mowgli_tea_enrich_pvals.pdf", bbox_inches="tight")
plt.close()


##########################################################################################
##################################### Annotate MOFA+ #####################################
##########################################################################################

# Compute neighbors UMAP embedding and Leiden clustering for MOFA+.
sc.pp.neighbors(mdata, n_neighbors=25, key_added="mofa_neighbors", use_rep="X_mofa")
sc.tl.umap(mdata, neighbors_key="mofa_neighbors")
sc.tl.leiden(
    mdata,
    resolution=0.2,
    key_added="leiden_mofa",
    neighbors_key="mofa_neighbors",
)

# Make the UMAP plot, colored by cluster.
ax = sc.pl.umap(
    mdata,
    color="leiden_mofa",
    alpha=0.7,
    legend_loc="on data",
    title="Leiden clustering of MOFA+ embedding",
    legend_fontoutline=2,
    frameon=False,
    show=False,
)
plt.savefig(fig_folder + "mofa_tea_umap.pdf", bbox_inches="tight")
plt.close()

# Make a dotplot with ADT values for each cluster.
mdata["adt"].obs["leiden_mofa"] = mdata.obs["leiden_mofa"]
mdata["adt"].var_names = mdata["adt"].var_names.str.replace("adt:", "")
axes = sc.pl.dotplot(
    mdata["adt"],
    adt_markers,
    groupby="leiden_mofa",
    title="MOFA+: counts for marker proteins in each cluster",
    categories_order=["9", "1", "8", "3", "5", "7", "6", "4", "0", "2"],
    mean_only_expressed=True,
    colorbar_title="Mean ADT count",
    size_title="Fraction of cells\nin cluster (%)",
    expression_cutoff=0.5,
    show=False,
)
axes["mainplot_ax"].set_ylabel("Cluster #")
axes["mainplot_ax"].set_xlabel("Marker proteins")
plt.savefig(fig_folder + "mofa_tea_leiden_adt.pdf", bbox_inches="tight")
plt.close()

# Make a dotplot with RNA values for each cluster.
rna_all_genes.obs["leiden_mofa"] = mdata.obs["leiden_mofa"]
axes = sc.pl.dotplot(
    rna_all_genes,
    rna_markers,
    categories_order=["9", "1", "8", "3", "5", "7", "6", "4", "0", "2"],
    groupby="leiden_mofa",
    title="Counts for marker genes in each cluster",
    mean_only_expressed=True,
    colorbar_title="Mean gene counts",
    size_title="Fraction of cells\nin cluster (%)",
    show=False,
)
axes["mainplot_ax"].set_ylabel("Cluster #")
axes["mainplot_ax"].set_xlabel("Marker genes")
plt.savefig(fig_folder + "mofa_tea_leiden_rna.pdf", bbox_inches="tight")
plt.close()

# Annotate the MOFA+ embedding.
mofa_cluster_names = {
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
}
codes = mdata.obs["leiden_mofa"].cat.codes
mdata.obs["annotation_mofa"] = [mofa_cluster_names[c] for c in codes]

# Make a UMAP plot of the MOFA+ embedding, colored by annotation.
ax = sc.pl.umap(
    mdata,
    color="annotation_mofa",
    alpha=0.7,
    legend_loc="on data",
    title="Annotated MOFA+ embedding",
    legend_fontoutline=2,
    frameon=False,
    legend_fontweight="normal",
    show=False,
)
plt.savefig(fig_folder + "mofa_tea_umap_annotated.pdf", bbox_inches="tight")
plt.close()

##########################################################################################
############################### Interpret MOFA's dimensions ##############################
##########################################################################################

# Make a dotplot of positive weights for MOFA+'s factors across clusters.
mofa_embedding = ad.AnnData(mdata.obsm["X_mofa"])
mofa_embedding.obs_names = mdata.obs_names
mofa_embedding.obs["annotation_mofa"] = mdata.obs["annotation_mofa"]
axes = sc.pl.dotplot(
    mofa_embedding,
    mofa_pos_factor_markers,
    groupby="annotation_mofa",
    categories_order=[
        "Erythroid cells",
        "B cells",
        "Monocytes",
        "NK cells",
        "MAIT T cells",
        "CD8 T cells",
        "CD4 T cells",
    ],
    expression_cutoff=0,
    mean_only_expressed=True,
    title="Positive weights of MOFA+'s factors",
    colorbar_title="Mean weight",
    size_title="% of cells with\nweights > 0 in cluster",
    cmap="Reds",
    show=False,
)
axes["mainplot_ax"].set_ylabel("Cluster")
axes["mainplot_ax"].set_xlabel("Factor #")
axes["mainplot_ax"].set_xticklabels(axes["mainplot_ax"].get_xticklabels(), rotation=0)
plt.savefig(fig_folder + "mofa_tea_factors_pos.pdf", bbox_inches="tight")
plt.close()


# Make a dotplot of negative weights for MOFA+'s factors across clusters.
mofa_embedding = ad.AnnData(-mdata.obsm["X_mofa"])
mofa_embedding.obs_names = mdata.obs_names
mofa_embedding.obs["annotation_mofa"] = mdata.obs["annotation_mofa"]
axes = sc.pl.dotplot(
    mofa_embedding,
    mofa_neg_factor_markers,
    groupby="annotation_mofa",
    categories_order=[
        "Erythroid cells",
        "B cells",
        "Monocytes",
        "NK cells",
        "MAIT T cells",
        "CD8 T cells",
        "CD4 T cells",
    ],
    expression_cutoff=0,
    mean_only_expressed=True,
    title="Negative weights of MOFA+'s factors",
    colorbar_title="Mean absolute weight",
    size_title="% of cells with\nweights < 0 in cluster",
    cmap="Blues",
    show=False,
)
axes["mainplot_ax"].set_ylabel("Cluster")
axes["mainplot_ax"].set_xlabel("Factor #")
axes["mainplot_ax"].set_xticklabels(axes["mainplot_ax"].get_xticklabels(), rotation=0)
plt.savefig(fig_folder + "mofa_tea_factors_neg.pdf", bbox_inches="tight")
plt.close()

# TODO: equal axis for bubble plots

# TODO: adt pos makers flat and adt neg markers flat
# TODO: standardize manually

# Make a matrixplot of positive ADT weights accross MOFA+'s factors.
adata = ad.AnnData(H_mofa["H_adt"].copy())
adata.X /= adata.X.std(axis=1, keepdims=True)
adata.X[adata.X < 0] = 0
adata.obs_names = mdata["adt"].var_names
adata.obs["adt"] = pd.Categorical(adata.obs_names)
adata = adata[adt_markers_flat, mofa_pos_factor_markers_flat]
axes = sc.pl.matrixplot(
    adata,
    mofa_pos_factor_markers,
    groupby="adt",
    cmap="Reds",
    categories_order=adt_markers_flat,
    title="Positive protein weights in MOFA+'s top factors",
    colorbar_title="ADT weight",
    # standard_scale="group",
    show=False,
)
axes["mainplot_ax"].set_ylabel("Proteins")
axes["mainplot_ax"].set_xlabel("Factor #")
plt.savefig(fig_folder + "mofa_tea_factors_pos_adt.pdf", bbox_inches="tight")
plt.close()

# Make a matrixplot of negative ADT weights accross MOFA+'s factors.
adata = ad.AnnData(-H_mofa["H_adt"].copy())
adata.X /= adata.X.std(axis=1, keepdims=True)
adata.X[adata.X < 0] = 0
adata.obs_names = mdata["adt"].var_names
adata.obs["adt"] = pd.Categorical(adata.obs_names)
adata = adata[adt_markers_flat, mofa_neg_factor_markers_flat]
axes = sc.pl.matrixplot(
    adata,
    mofa_neg_factor_markers,
    groupby="adt",
    cmap="Blues",
    categories_order=adt_markers_flat,
    title="Negative protein weights in MOFA+'s top factors",
    colorbar_title="Absolute ADT weight",
    # standard_scale="group",
    show=False,
)
axes["mainplot_ax"].set_ylabel("Proteins")
axes["mainplot_ax"].set_xlabel("Factor #")
plt.savefig(fig_folder + "mofa_tea_factors_neg_adt.pdf", bbox_inches="tight")
plt.close()

# Make a matrixplot of positive RNA weights accross MOFA+'s factors.
adata = ad.AnnData(H_mofa["H_rna"].copy())
adata.X /= adata.X.std(axis=1, keepdims=True)
adata.X[adata.X < 0] = 0
adata.obs_names = mdata["rna"].var_names.str.replace("rna:", "")
adata.obs["rna"] = pd.Categorical(adata.obs_names)
genes = [g for g in rna_markers_flat if g in adata.obs_names]
adata = adata[genes, mofa_pos_factor_markers_flat]
sc.pl.matrixplot(
    adata,
    mofa_pos_factor_markers,
    groupby="rna",
    cmap="Reds",
    categories_order=genes,
    title="Positive gene weights in MOFA+'s top factors",
    colorbar_title="Gene weight",
    # standard_scale="group",
    show=False,
)
plt.savefig(fig_folder + "mofa_tea_factors_pos_rna.pdf", bbox_inches="tight")
plt.close()

# Make a matrixplot of negative RNA weights accross MOFA+'s factors.
adata = ad.AnnData(-H_mofa["H_rna"].copy())
adata.X /= adata.X.std(axis=1, keepdims=True)
adata.X[adata.X < 0] = 0
adata.obs_names = mdata["rna"].var_names.str.replace("rna:", "")
adata.obs["rna"] = pd.Categorical(adata.obs_names)
genes = [g for g in rna_markers_flat if g in adata.obs_names]
adata = adata[genes, mofa_neg_factor_markers_flat]
sc.pl.matrixplot(
    adata,
    mofa_neg_factor_markers,
    groupby="rna",
    cmap="Blues",
    categories_order=genes,
    title="Negative gene weights in MOFA+'s top factors",
    show=False,
    colorbar_title="Absolute gene weight",
    # standard_scale="group",
)
plt.savefig(fig_folder + "mofa_tea_factors_neg_rna.pdf", bbox_inches="tight")
plt.close()

mofa_pvals_top = ad.AnnData(np.zeros((len(gene_sets), mofa_embedding.n_vars)))
mofa_pvals_top.obs_names = gene_sets.keys()
mofa_pvals_top.var_names = mofa_embedding.var_names

for i, x in enumerate(gene_sets):
    idx = enr["native"] == gene_sets[x]
    idx &= enr["query"].str.startswith("top_mofa")
    mofa_pvals_top.X[i, enr.loc[idx, "dim"]] = enr.loc[idx, "minlogp"]

mofa_pvals_top.obs["gene_set"] = pd.Categorical(mofa_pvals_top.obs_names)
axes = sc.pl.matrixplot(
    mofa_pvals_top,
    mofa_pos_factor_markers,
    groupby="gene_set",
    categories_order=gene_sets.keys(),
    cmap="Reds",
    title="MOFA+: Cell-type gene set enrichment for the top genes in factors",
    colorbar_title=r"$-\log_{10}$(p value)",
    show=False,
)
axes["mainplot_ax"].set_xlabel("Factor #")
axes["mainplot_ax"].set_ylabel("Gene set")
plt.savefig(fig_folder + "mofa_tea_enrich_pvals_pos.pdf", bbox_inches="tight")
plt.close()

mofa_pvals_bottom = ad.AnnData(np.zeros((len(gene_sets), mofa_embedding.n_vars)))
mofa_pvals_bottom.obs_names = gene_sets.keys()
mofa_pvals_bottom.var_names = mofa_embedding.var_names

for i, x in enumerate(gene_sets):
    idx = enr["native"] == gene_sets[x]
    idx &= enr["query"].str.startswith("bottom_mofa")
    mofa_pvals_bottom.X[i, enr.loc[idx, "dim"]] = enr.loc[idx, "minlogp"]

mofa_pvals_bottom.obs["gene_set"] = pd.Categorical(mofa_pvals_bottom.obs_names)
axes = sc.pl.matrixplot(
    mofa_pvals_bottom,
    mofa_neg_factor_markers,
    groupby="gene_set",
    categories_order=gene_sets.keys(),
    cmap="Blues",
    title="MOFA+: Cell-type gene set enrichment for the bottom genes in factors",
    colorbar_title=r"$-\log_{10}$(p value)",
    show=False,
)
axes["mainplot_ax"].set_xlabel("Factor #")
axes["mainplot_ax"].set_ylabel("Gene set")
plt.savefig(fig_folder + "mofa_tea_enrich_pvals_neg.pdf", bbox_inches="tight")
plt.close()

##########################################################################################
###################################### Bubble plots ######################################
##########################################################################################

# TODO: same sizes everywhere. Some sort of score.
# TODO: Ratio
# TODO: Numbers of the factors

mowgli_embedding = ad.AnnData(mdata.obsm["X_mowgli"])
mowgli_embedding.obs_names = mdata.obs_names
mowgli_embedding.obs["annotation_mowgli"] = mdata.obs["annotation_mowgli"]
fig, axes = plt.subplots(2, 7, figsize=(25, 6), constrained_layout=True)

for i, cluster in enumerate(
    [
        "Erythroid cells",
        "B cells",
        "Monocytes",
        "NK cells",
        "MAIT T cells",
        "CD8 T cells",
        "CD4 T cells",
    ]
):

    df = []

    for factor in mowgli_embedding.var_names:
        idx_cluster = mowgli_embedding.obs["annotation_mowgli"] == cluster
        mean_cluster = float(mowgli_embedding[idx_cluster, factor].X.ravel().mean())
        mean_other = float(mowgli_embedding[~idx_cluster, factor].X.ravel().mean())
        prop_expressed = float(
            np.mean(mowgli_embedding[idx_cluster, factor].X.ravel() > 1e-3)
        )
        df.append(
            {
                "factor": factor,
                "mean_cluster": mean_cluster,
                "mean_other": mean_other,
                "prop_expressed": prop_expressed,
            }
        )
    df = pd.DataFrame(df)

    df = df[df["prop_expressed"] > 0]

    max_val = df[["mean_other", "mean_cluster"]].to_numpy().max()
    padding = 0.1 * max_val
    axes[0, i].plot(
        [0, max_val], [0, max_val], linewidth=2, linestyle="--", color="grey"
    )
    sns.scatterplot(
        data=df,
        x="mean_other",
        y="mean_cluster",
        size="prop_expressed",
        ax=axes[0, i],
        sizes=(500 * df["prop_expressed"].min(), 500 * df["prop_expressed"].max()),
        alpha=0.7,
        color="purple",
    )

    for j in range(len(df)):
        if (df["mean_cluster"].iloc[j] - df["mean_other"].iloc[j]) / df[
            "mean_cluster"
        ].max() > 0.25:
            axes[0, i].text(
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

    j = np.argmax((df["mean_cluster"] - df["mean_other"]) / df["mean_cluster"].max())
    score = (df["mean_cluster"].iloc[j] - df["mean_other"].iloc[j]) / df[
        "mean_cluster"
    ].max()
    axes[0, i].text(
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

    if i < 6:
        axes[0, i].get_legend().remove()
    else:
        axes[0, i].legend(
            loc="upper left",
            bbox_to_anchor=(1.05, 1),
            frameon=False,
            title="Prop. of cells in cluster\nwith weights > 1e-3",
        )
        for lh in axes[0, i].get_legend().legendHandles:
            lh.set_alpha(0.7)

    axes[0, i].set_title(cluster)
    axes[0, i].set_xlabel(None)
    if i > 0:
        axes[0, i].set_ylabel(None)
    else:
        axes[0, i].set_ylabel("Mean weight in cell type")
    axes[0, i].set_xlim(-padding, max_val + padding)
    axes[0, i].set_ylim(-padding, max_val + padding)

    axes[0, i].spines.right.set_visible(False)
    axes[0, i].spines.top.set_visible(False)
    axes[0, i].yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    axes[0, i].xaxis.set_major_formatter(FormatStrFormatter("%.2f"))

mofa_embedding_pos = ad.AnnData(mdata.obsm["X_mofa"])
mofa_embedding_pos.obs_names = mdata.obs_names
mofa_embedding_pos.obs["annotation_mofa"] = mdata.obs["annotation_mofa"]
mofa_embedding_pos.X[mofa_embedding_pos.X < 0] = 0

mofa_embedding_neg = ad.AnnData(-mdata.obsm["X_mofa"])
mofa_embedding_neg.obs_names = mdata.obs_names
mofa_embedding_neg.obs["annotation_mofa"] = mdata.obs["annotation_mofa"]
mofa_embedding_neg.X[mofa_embedding_neg.X < 0] = 0

for i, cluster in enumerate(
    [
        "Erythroid cells",
        "B cells",
        "Monocytes",
        "NK cells",
        "MAIT T cells",
        "CD8 T cells",
        "CD4 T cells",
    ]
):

    df = []

    for factor in mofa_embedding_pos.var_names:
        idx_cluster = mofa_embedding_pos.obs["annotation_mofa"] == cluster
        mean_cluster = float(mofa_embedding_pos[idx_cluster, factor].X.ravel().mean())
        mean_other = float(mofa_embedding_pos[~idx_cluster, factor].X.ravel().mean())
        prop_expressed = float(
            np.mean(mofa_embedding_pos[idx_cluster, factor].X.ravel() > 1e-3)
        )
        df.append(
            {
                "factor": factor,
                "Sign of weights": "Positive",
                "mean_cluster": mean_cluster,
                "mean_other": mean_other,
                "Proportion": prop_expressed,
            }
        )

    for factor in mofa_embedding_neg.var_names:
        idx_cluster = mofa_embedding_neg.obs["annotation_mofa"] == cluster
        mean_cluster = float(mofa_embedding_neg[idx_cluster, factor].X.ravel().mean())
        mean_other = float(mofa_embedding_neg[~idx_cluster, factor].X.ravel().mean())
        prop_expressed = float(
            np.mean(mofa_embedding_neg[idx_cluster, factor].X.ravel() > 1e-3)
        )
        df.append(
            {
                "factor": factor,
                "Sign of weights": "Negative",
                "mean_cluster": mean_cluster,
                "mean_other": mean_other,
                "Proportion": prop_expressed,
            }
        )

    df = pd.DataFrame(df)

    max_val = df[["mean_other", "mean_cluster"]].to_numpy().max()
    padding = 0.1 * max_val
    axes[1, i].plot(
        [0, max_val], [0, max_val], linewidth=2, linestyle="--", color="grey"
    )

    sns.scatterplot(
        data=df,
        x="mean_other",
        y="mean_cluster",
        hue="Sign of weights",
        palette=["tab:blue", "tab:red"],
        size="Proportion",
        ax=axes[1, i],
        sizes=(
            500 * df["Proportion"].min(),
            500 * df["Proportion"].max(),
        ),
        alpha=0.7,
    )

    for j in range(len(df)):
        if (df["mean_cluster"].iloc[j] - df["mean_other"].iloc[j]) / df[
            "mean_cluster"
        ].max() > 0.25:
            axes[1, i].text(
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

    j = np.argmax((df["mean_cluster"] - df["mean_other"]) / df["mean_cluster"].max())
    score = (df["mean_cluster"].iloc[j] - df["mean_other"].iloc[j]) / df[
        "mean_cluster"
    ].max()
    axes[1, i].text(
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

    if i > 0:
        axes[1, i].set_ylabel(None)
    else:
        axes[1, i].set_ylabel("Mean weight in cell type")

    if i < 6:
        axes[1, i].get_legend().remove()
    else:
        axes[1, i].legend(
            loc="upper left",
            bbox_to_anchor=(1.05, 1),
            frameon=False,
            title="Prop. of cells in cluster with\npositive / negative weights",
        )
        for lh in axes[1, i].get_legend().legendHandles:
            lh.set_alpha(0.7)

    axes[1, i].spines.right.set_visible(False)
    axes[1, i].spines.top.set_visible(False)
    axes[1, i].yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    axes[1, i].xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    axes[1, i].set_xlabel("Mean weight in other cell types")
    axes[1, i].set_xlim(-padding, max_val + padding)
    axes[1, i].set_ylim(-padding, max_val + padding)

plt.savefig(fig_folder + "tea_factors_bubble.pdf", bbox_inches="tight")
plt.close()

##########################################################################################
############################### Focus: Memory B cell factor ##############################
##########################################################################################

fig, axes = plt.subplots(1, 3, figsize=(10, 3), constrained_layout=True)
sc.pp.neighbors(mowgli_embedding, n_neighbors=25)
sc.tl.umap(mowgli_embedding)
sc.pl.embedding(
    mowgli_embedding,
    color="44",
    alpha=0.7,
    title="Weight of factor 44",
    basis="X_umap",
    frameon=False,
    show=False,
    ax=axes[0],
)
df = pd.DataFrame(
    H_mowgli["H_adt"],
    index=mdata["adt"].var_names,
    columns=mowgli_embedding.var_names,
)

sns.barplot(
    data=df,
    y="44",
    x=df.index,
    ax=axes[1],
    order=df.index[df["44"].argsort()[::-1][:15]],
    color="purple",
)
axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=90)

sns.barplot(
    data=enr[
        (enr["method"] == "mowgli")
        & (enr["dim"] == 44)
        & (enr["native"].str.contains("GO"))
    ][:10],
    x="minlogp",
    y="native",
    ax=axes[2],
)
# Shorten the yticklabels to 15 characters.
axes[2].set_yticklabels([label.get_text()[:15] for label in axes[2].get_yticklabels()])
# axes[2].set_yticks(axes[2].get_yticks(), [label.get_text()[:15] for label in axes[2].get_yticklabels()])

plt.savefig(fig_folder + "mowgli_tea_umap_44.pdf", bbox_inches="tight")
plt.close()
