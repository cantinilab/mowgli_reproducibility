######################################### Imports ########################################

import anndata as ad
import matplotlib.pyplot as plt
import mofax
import muon as mu
import numpy as np
import os
import pandas as pd
import plotly.graph_objects as go
import scanpy as sc
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scanpy.plotting import palettes
from urllib.error import HTTPError

#################################### Define some paths ###################################

# Define the data, enrichment, and figures folders.
data_folder = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/"
enrich_folder = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/enrich/"
fig_folder = os.path.join(
    "/users/csb/huizing/Documents/PhD/Code/",
    "mowgli_reproducibility/visualize/figures/",
)

# Define the folder containing Mowgli's embeddings and loadings.
w_folder = "/users/csb/huizing/Documents/PhD/Code/Mowgli/local_analysis/from_jz/w/"
h_folder = "/users/csb/huizing/Documents/PhD/Code/Mowgli/local_analysis/from_jz/h/"

#################################### Load the dataset ####################################

# Load the preprocessed data.
mdata = mu.read_h5mu(os.path.join(data_folder, "TEA/tea_preprocessed.h5mu.gz"))
mdata.uns = {}  # This is needed somehow.
mdata.obs[mdata["adt"].var_names] = mdata["adt"].X

################################## Load the Mowgli model #################################

# Load the Mowgli model.
mdata.obsm["X_mowgli"] = np.load(
    os.path.join(
        w_folder, "tea_mowgli_cosine_50_0_05_rna_0_01_atac_0_1_adt_0_01_0_001.npy"
    ),
    allow_pickle=True,
).item()["W"]

# Load Mowgli's loadings.
H_mowgli = np.load(
    os.path.join(
        h_folder, "tea_mowgli_cosine_50_0_05_rna_0_01_atac_0_1_adt_0_01_0_001.npy"
    ),
    allow_pickle=True,
).item()

# Make an AnnData object for adt.
H_adt = ad.AnnData(H_mowgli["H_adt"])
H_adt.X /= H_adt.X.sum(axis=1, keepdims=True)
H_adt.obs_names = mdata["adt"].var_names.str.replace("adt:", "")

# Make an AnnData object for rna.
H_rna = ad.AnnData(H_mowgli["H_rna"])
H_rna.X /= H_rna.X.sum(axis=1, keepdims=True)
H_rna.obs_names = mdata["rna"].var_names.str.replace("rna:", "")

# Add mowgli dimensions.
mdata.obs[[f"mowgli:{i}" for i in range(50)]] = mdata.obsm["X_mowgli"]

################################### Load the MOFA model ##################################

# Load the MOFA model.
mofa_model = mofax.mofa_model(os.path.join(data_folder, "TEA/tea_mofa_15.hdf5"))
mdata.obsm["X_mofa"] = mofa_model.get_factors()

# Add mofa dimensions.
mdata.obs[[f"mofa:{i}" for i in range(15)]] = mdata.obsm["X_mofa"]

################################# Load the azimuth annotation ############################

# Load the azimuth annotation.
azimuth_pred = pd.read_csv(
    os.path.join(data_folder, "TEA/azimuth_pred.tsv"),
    sep="\t",
    index_col=0,
)

mdata.obs[azimuth_pred.columns] = azimuth_pred.loc[mdata.obs_names]

############################## Load the gene enrichment results ##########################

# Load the gprofiler results, rename columns and make dimensions integers.
enr = pd.read_csv(os.path.join(enrich_folder, "enrichment.csv"), index_col=0)
enr["dim"] = enr["query"].str.split(" ", expand=True)[1].astype(int)

# Add some more columns.
enr["name"] = enr["native"]
enr = enr[enr["p_value"] < 0.05]
enr["minlogp"] = -np.log10(enr["p_value"])

############################ Load the motif enrichment results ###########################

# Define the motifs.
all_motifs = pd.DataFrame({})
for i in range(50):
    motifs = pd.read_csv(
        os.path.join(enrich_folder, f"top_motifs_mowgli/motifs_{i}.csv"),
        index_col=0,
    )
    motifs = motifs[motifs["p.adjust"] < 0.05]
    motifs = motifs[motifs["motif"].str.startswith("MA")]
    motifs["dim"] = i
    motifs["method"] = "mowgli"
    all_motifs = pd.concat((all_motifs, motifs))

all_motifs["minlogp"] = -np.log10(all_motifs["p.adjust"])

####################### Compute UMAP embeddings and Leiden clusters ######################

# Compute the UMAP embedding and Leiden clusters for MOFA.
sc.pp.neighbors(mdata, n_neighbors=25, key_added="mofa_neighbors", use_rep="X_mofa")
sc.tl.leiden(mdata, 0.2, key_added="mofa_leiden", neighbors_key="mofa_neighbors")
sc.tl.umap(mdata, neighbors_key="mofa_neighbors")
mdata.obsm["X_umap_mofa"] = mdata.obsm["X_umap"]

# Compute the UMAP embedding and Leiden clusters for Mowgli.
sc.pp.neighbors(mdata, n_neighbors=25, key_added="mowgli_neighbors", use_rep="X_mowgli")
sc.tl.leiden(mdata, 0.2, key_added="mowgli_leiden", neighbors_key="mowgli_neighbors")
sc.tl.umap(mdata, neighbors_key="mowgli_neighbors")
mdata.obsm["X_umap_mowgli"] = mdata.obsm["X_umap"]

################################## Name MOFA's clusters ##################################

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
codes = mdata.obs["mofa_leiden"].cat.codes
mdata.obs["annotation_mofa"] = [mofa_cluster_names[c] for c in codes]

################################## Name Mowgli's clusters ################################

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
codes = mdata.obs["mowgli_leiden"].cat.codes
mdata.obs["annotation_mowgli"] = [mowgli_cluster_names[c] for c in codes]

######################### Plot the UMAPs with manual annotation ##########################

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

sc.pl.embedding(
    mdata,
    basis="X_umap_mofa",
    color="annotation_mofa",
    title="MOFA-based annotation",
    legend_loc="on data",
    legend_fontoutline=2,
    alpha=0.7,
    frameon=False,
    ax=axes[0],
    show=False,
)

sc.pl.embedding(
    mdata,
    basis="X_umap_mowgli",
    color="annotation_mowgli",
    title="Mowgli-based annotation",
    legend_loc="on data",
    legend_fontoutline=2,
    alpha=0.7,
    frameon=False,
    ax=axes[1],
    show=False,
)
plt.savefig(os.path.join(fig_folder, "annotation_manual.pdf"), bbox_inches="tight")

######################## Plot the UMAPs with azimuth annotation ##########################

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

sc.pl.embedding(
    mdata,
    basis="X_umap_mofa",
    color="predicted.celltype.l3",
    title="Azimuth annotation",
    legend_loc="none",
    alpha=0.7,
    frameon=False,
    ax=axes[0],
    show=False,
    palette=palettes.default_102,
)
sc.pl.embedding(
    mdata,
    basis="X_umap_mowgli",
    color="predicted.celltype.l3",
    title="Azimuth annotation",
    frameon=False,
    alpha=0.7,
    ax=axes[1],
    show=False,
    palette=palettes.default_102,
)
plt.savefig(os.path.join(fig_folder, "annotation_azimuth.pdf"), bbox_inches="tight")

##################################### MOFA embedding #####################################

sc.pl.embedding(
    mdata,
    basis="X_umap_mofa",
    alpha=0.7,
    color=[f"mofa:{i}" for i in range(15)],
    legend_fontoutline=2,
    frameon=False,
    color_map="RdBu",
    vmax=lambda x: np.max(np.abs(x)),
    vmin=lambda x: -np.max(np.abs(x)),
    show=False,
)
plt.savefig(os.path.join(fig_folder, "mofa.pdf"), bbox_inches="tight")

#################################### Mowgli embedding ####################################

sc.pl.embedding(
    mdata,
    basis="X_umap_mowgli",
    alpha=0.7,
    color=[f"mowgli:{i}" for i in range(50)],
    legend_fontoutline=2,
    frameon=False,
    color_map="Purples",
    show=False,
)
plt.savefig(os.path.join(fig_folder, "mowgli.pdf"), bbox_inches="tight")

########################### Mowgli embedding for the proteins ############################

sc.pl.embedding(
    mdata,
    basis="X_umap_mowgli",
    alpha=0.7,
    color=mdata["adt"].var_names.to_list(),
    legend_fontoutline=2,
    frameon=False,
    show=False,
    save=False,
)
plt.savefig(os.path.join(fig_folder, "mowgli_adt.pdf"), bbox_inches="tight")

############################ MOFA embedding for the proteins #############################

sc.pl.embedding(
    mdata,
    basis="X_umap_mofa",
    alpha=0.7,
    color=mdata["adt"].var_names.to_list(),
    legend_fontoutline=2,
    frameon=False,
    show=False,
    save=False,
)
plt.savefig(os.path.join(fig_folder, "mofa_adt.pdf"), bbox_inches="tight")

####################### Sankey diagram of the cell type annotations ######################

annot_azimuth = mdata.obs["predicted.celltype.l1"]
annot_mofa = mdata.obs["annotation_mofa"]
annot_mowgli = mdata.obs["annotation_mowgli"]

categories_mowgli = np.sort(annot_mowgli.cat.categories).tolist()
categories_azimuth = np.sort(annot_azimuth.cat.categories).tolist()
categories_mofa = np.sort(annot_mofa.cat.categories).tolist()

categories = categories_mowgli + categories_azimuth + categories_mofa
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

colors = ["#002626", "#0E4749", "#E6E075", "#F07C42", "#C2A370", "#142A58"]
colors += ["#C87E7E", "#002626", "#0E4749", "#E6E075", "#F36C6C", "#142A58"]
colors += ["#C87E7E", "#F07C42", "#C2A370", "#002626", "#0E4749", "#E6E075"]
colors += ["#F07C42", "#C2A370", "#142A58", "#C87E7E"]

line = {"color": "black", "width": 0.5}
node = {"pad": 10, "thickness": 20, "line": line, "label": categories, "color": colors}
link = {"source": source, "target": target, "value": value}

fig = go.Figure(data=[go.Sankey(node=node, link=link)])
fig.update_layout(font_size=10, width=500, height=500)
fig.write_image(os.path.join(fig_folder, "sankey.pdf"))
# TODO: how to save this?

############################### Translate proteins to genes ##############################

adt_to_rna = {"CD10": "MME", "CD11b": "ITGAM", "CD11c": "ITGAX", "CD123": "IL3RA"}
adt_to_rna = {**adt_to_rna, "CD127": "IL7R", "CD14": "CD14", "CD141": "THBD"}
adt_to_rna = {**adt_to_rna, "CD172a": "SIRPA", "CD185": "CXCR5", "CD19": "CD19"}
adt_to_rna = {**adt_to_rna, "CD197": "CCR7", "CD21": "CR2", "CD24": "CD24"}
adt_to_rna = {**adt_to_rna, "CD269": "TNFRSF17", "CD27": "CD27", "CD278": "ICOS"}
adt_to_rna = {**adt_to_rna, "CD3": "CD3G", "CD304": "NRP1", "CD319": "SLAMF7"}
adt_to_rna = {**adt_to_rna, "CD39": "ENTPD1", "CD4": "CD4", "CD40": "CD40"}
adt_to_rna = {**adt_to_rna, "CD45RO": "PTPRC", "CD56": "NCAM1", "CD66b": "CEACAM8"}
adt_to_rna = {**adt_to_rna, "CD80": "CD80", "CD86": "CD86", "CD8a": "CD8A"}
adt_to_rna = {**adt_to_rna, "HLA-DR": "HLA-DRA", "IgD": "IGHD", "IgM": "IGHM"}
adt_to_rna = {**adt_to_rna, "CD16": "FCGR3A", "CD192": "CCR2", "CD25": "IL2RA"}
adt_to_rna = {**adt_to_rna, "CD279": "PDCD1", "CD38": "CD38", "CD45RA": "PTPRC"}
adt_to_rna = {**adt_to_rna, "CD71": "TFRC", "CD95": "FAS", "KLRG1": "KLRG1"}
# "FceRI": "FceRI"
# "IgG1-K-Isotype-Control": "IgG1-K-Isotype-Control"
# "TCR-Va24-Ja18": "TCR-Va24-Ja18"
# "TCR-Va7.2": "TCR-Va7.2"
# "TCR-a/b": "TCR-a/b"
# "TCR-g/d": "TCR-g/d"

########################## Are top genes targets of top motifs? ##########################


def get_tf_gene_matches(dim: int, pathway_file: str, n_genes=20, n_tf=20):

    pathways = pd.read_csv(
        os.path.join(data_folder, "pathways/", pathway_file),
        sep="\t",
        header=None,
        names=["TF", "target", "weight"],
    ).sort_values("weight", ascending=False)

    idx = np.argsort(H_rna.X[:, dim])[::-1][:n_genes]
    top_rna = H_rna.obs_names[idx].to_list()

    top_motifs = all_motifs[all_motifs["dim"] == dim].sort_values("pvalue")
    top_motifs = top_motifs.head(n_tf)["motif.name"].to_list()

    list_motifs = []
    for motif in top_motifs:
        idx = pathways["TF"].apply(lambda g: g in motif.split("::"))
        links = pathways.loc[idx, "target"].to_list()
        if len(np.intersect1d(links, top_rna)) > 0:
            list_motifs.append(motif)

    return list_motifs


###################### Are top genes or proteins cell-type specific? #####################


def add_stars(gene_names, cell_types, adt=False):

    updated_gene_names = []
    link = "https://www.proteinatlas.org/search/"
    for gene_name in gene_names:

        try:
            if adt and gene_name in adt_to_rna:
                search_name = adt_to_rna[gene_name]
            else:
                search_name = gene_name
            res = pd.read_csv(link + f"{search_name}?format=tsv", sep="\t")
        except HTTPError:
            res = pd.DataFrame({})

        specific = {}
        if len(res) > 0 and not pd.isna(res.loc[0, "RNA blood cell specific nTPM"]):
            for x in res.loc[0, "RNA blood cell specific nTPM"].split(";"):
                name, score = x.split(": ")
                specific[name] = float(score)

        new_gene_name = gene_name
        for cell_type in cell_types:
            if cell_type in specific:
                new_gene_name = f"* {gene_name}"
                break

        updated_gene_names.append(new_gene_name)

    assert len(gene_names) == len(updated_gene_names)
    return updated_gene_names


############################# Make the big plot with 4 zooms #############################


fig, ax = plt.subplots(1, 1, figsize=(7, 4))
sc.pl.umap(
    mdata,
    alpha=0.7,
    legend_loc=None,
    legend_fontoutline=2,
    frameon=False,
    show=False,
    ax=ax,
)

idx_zoom_1 = mdata.obs["annotation_mowgli"].str.contains("NK")
axins_1 = inset_axes(
    ax,
    "40%",
    "60%",
    borderpad=3,
    loc="upper left",
    bbox_to_anchor=(-0.1, -0.6, 1, 1),
    bbox_transform=ax.transAxes,
)
sc.pl.umap(
    mdata[idx_zoom_1],
    alpha=0.7,
    color=["mowgli:2"],
    legend_loc=None,
    legend_fontoutline=2,
    vmax=mdata.obs["mowgli:2"].max(),
    show=False,
    ax=axins_1,
)
axins_1.set_title("CD56-dim Natural Killer Cells", fontsize=10)
fig.axes[-1].remove()

idx_zoom_2 = mdata.obs["annotation_mowgli"].str.contains("B")
axins_2 = inset_axes(
    ax,
    "40%",
    "60%",
    borderpad=3,
    loc="upper left",
    bbox_to_anchor=(0.75, -0.6, 1, 1),
    bbox_transform=ax.transAxes,
)
sc.pl.umap(
    mdata[idx_zoom_2],
    alpha=0.7,
    color=["mowgli:33"],
    legend_loc=None,
    legend_fontoutline=2,
    vmax=mdata.obs["mowgli:33"].max(),
    show=False,
    ax=axins_2,
)
axins_2.set_title("Naive B cells", fontsize=10)
fig.axes[-1].remove()

idx_zoom_3 = mdata.obs["annotation_mowgli"].str.contains("B")
axins_3 = inset_axes(
    ax,
    "40%",
    "60%",
    borderpad=3,
    loc="upper left",
    bbox_to_anchor=(0.75, 0.2, 1, 1),
    bbox_transform=ax.transAxes,
)
sc.pl.umap(
    mdata[idx_zoom_3],
    alpha=0.7,
    color=["mowgli:44"],
    legend_loc=None,
    legend_fontoutline=2,
    vmax=mdata.obs["mowgli:44"].max(),
    show=False,
    ax=axins_3,
)
axins_3.set_title("Memory B cells", fontsize=10)
fig.axes[-1].remove()

idx_zoom_4 = mdata.obs["annotation_mowgli"].str.contains("CD8")
axins_4 = inset_axes(
    ax,
    "40%",
    "60%",
    borderpad=3,
    loc="upper left",
    bbox_to_anchor=(-0.1, 0.2, 1, 1),
    bbox_transform=ax.transAxes,
)
sc.pl.umap(
    mdata[idx_zoom_4],
    alpha=0.7,
    color=["mowgli:49"],
    legend_loc=None,
    legend_fontoutline=2,
    vmax=mdata.obs["mowgli:49"].max(),
    show=False,
    ax=axins_4,
)
axins_4.set_title("Effector Memory CD8 T cells", fontsize=10)
fig.axes[-1].remove()

for axins, idx_zoom in zip(
    [axins_1, axins_2, axins_3, axins_4],
    [idx_zoom_1, idx_zoom_2, idx_zoom_3, idx_zoom_4],
):
    obsm = mdata[idx_zoom].obsm["X_umap"]
    x1, x2 = np.percentile(obsm[:, 0], 1), np.percentile(obsm[:, 0], 99)
    y1, y2 = np.percentile(obsm[:, 1], 1), np.percentile(obsm[:, 1], 99)
    padding_x, padding_y = (x2 - x1) * 0.3, (y2 - y1) * 0.3
    axins.set_xlim(x1 - padding_x, x2 + padding_x)
    axins.set_ylim(y1 - padding_y, y2 + padding_y)
    axins.set_xlabel(None)
    axins.set_ylabel(None)
    axins.set_facecolor((1, 1, 1, 0.75))
    rect, lines = ax.indicate_inset_zoom(axins)
    for line in lines:
        line.set_visible(True)

plt.savefig(os.path.join(fig_folder, "multi_umap.pdf"), bbox_inches="tight")

########################### Make a plot with dictionary weights ##########################


def plot_weights(dim, H, manual_stars, star_cell_types, n_features, ax, adt=False):
    idx = np.argsort(H.X[:, dim])[::-1][:n_features]
    names = H.obs_names[idx].to_list()
    names = add_stars(names, star_cell_types, adt=adt)
    for name in manual_stars:
        if name in names:
            names[names.index(name)] = f"* {name}"
    sns.barplot(y=names, x=H.X[idx, dim], ax=ax, palette="Blues_r")
    for direction in ["top", "right", "left"]:
        ax.spines[direction].set_visible(False)


################################ Make a plot with gene sets ##############################


def plot_gene_sets(enr, dim, contains, contains_not, n_gs, gs_stars, ax):
    idx_immune = enr["source"] == "immune"
    for word in contains:
        idx_immune &= enr["name"].str.split("VS", expand=True)[1].str.contains(word)
    for word in contains_not:
        idx_immune &= ~(enr["name"].str.contains(word))
    idx_azimuth = enr["source"] == "azimuth"
    idx_bp = enr["source"] == "go_bp"
    idx_cc = enr["source"] == "go_cc"
    idx_mf = enr["source"] == "go_mf"
    idx = (idx_immune | idx_azimuth | idx_bp | idx_cc | idx_mf) & (enr["dim"] == dim)

    df = enr[idx].sort_values("minlogp", ascending=False)

    idx = df["source"] == "immune"
    df.loc[idx, "name"] = (
        df.loc[idx, "name"]
        .apply(lambda x: " ".join(x.split("_")[1:]))
        .str.replace(" DN", "  down")
        .str.replace(" UP", "  up")
        .str.replace("VS", " vs ")
    )

    idx, pattern = (df["source"] == "azimuth"), " CL[0-9]+"
    df.loc[idx, "name"] = df.loc[idx, "name"].str.replace(pattern, "", regex=True)

    for go_source in ["go_bp", "go_cc", "go_mf"]:
        idx, pattern = (df["source"] == go_source), " \(GO:[0-9]+\)"
        df.loc[idx, "name"] = df.loc[idx, "name"].str.replace(pattern, "", regex=True)

    df = df.drop_duplicates("name").head(n_gs)

    names = df["name"].to_list()
    for i in gs_stars:
        names[i] = f"* {names[i]}"
    df["name"] = names

    ax.hlines(y=df["name"], xmin=0, xmax=df["minlogp"], color="lightgray")
    sns.scatterplot(
        data=df,
        x="minlogp",
        y="name",
        ax=ax,
        hue="source",
        hue_order=["azimuth", "immune", "go_bp", "go_cc", "go_mf"],
        palette="tab10",
        alpha=0.85,
        estimator=max,
        ci=None,
        s=100,
        zorder=3,
    )
    ax.set_xlabel(r"$-\log_{10}~p$")
    ax.set_ylabel(None)
    ax.set_title("Gene set enrichment")
    for direction in ["top", "right", "left"]:
        ax.spines[direction].set_visible(False)


################################# Make a plot with motifs ################################


def plot_motifs(all_motifs, dim, n_motifs, pathway_files, ax):
    df = all_motifs[all_motifs["dim"] == dim].sort_values("pvalue", ascending=True)
    df = df.drop_duplicates("motif.name").head(n_motifs)

    for pathway_file in pathway_files:
        list_motifs = get_tf_gene_matches(dim, pathway_file)
        for motif in list_motifs:
            df.loc[df["motif.name"] == motif, "motif.name"] = f"* {motif}"

    ax.hlines(y=df["motif.name"], xmin=0, xmax=df["minlogp"], color="lightgray")
    x, y = "minlogp", "motif.name"
    sns.scatterplot(data=df, y=y, x=x, ax=ax, palette="Blues_r", s=100, zorder=3)
    ax.set_xlabel(r"$-\log_{10}~p$")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_title("Motif enrichment")
    ax.set_ylabel(None)


########################## Assemble the previous plots into one ##########################


def compact_plot(
    dim,
    contains,
    contains_not,
    star_cell_types,
    manual_stars,
    pathway_files,
    gs_stars,
    axes,
    n_rna=15,
    n_adt=15,
    n_gs=20,
    n_motifs=20,
):

    plot_gene_sets(enr, dim, contains, contains_not, n_gs, gs_stars, axes[0])
    plot_motifs(all_motifs, dim, n_motifs, pathway_files, axes[1])

    plot_weights(dim, H_rna, manual_stars, star_cell_types, n_rna, axes[2], adt=False)
    axes[2].set_title("RNA weights")
    axes[2].set_xlabel("RNA weights")

    plot_weights(dim, H_adt, manual_stars, star_cell_types, n_adt, axes[3], adt=True)
    axes[3].set_title("ADT weights")
    axes[3].set_xlabel("ADT weights")

    for ax in fig.axes:
        for tick in ax.get_yticklabels():
            tick.set_fontname("Roboto Condensed")


################## Plot the weights en enrichments for different factors #################

fig, axes = plt.subplots(4, 4, figsize=(14, 10), constrained_layout=True)

compact_plot(
    dim=2,
    contains=["NK"],
    contains_not=["TCELL", "BCELL", "MONO", "DC"],
    star_cell_types=["NK-cell"],
    manual_stars=["CD11c", "KLRG1", "CD16"],
    gs_stars=[0, 1, 2, 4, 5, 6, 7, 8, 13, 14],  # , 15, 18, 19]
    pathway_files=["natural_killer_cells.txt.gz"],
    axes=axes[0],
)

compact_plot(
    dim=33,
    contains=["BCELL"],
    contains_not=["TCELL", "DC", "MONO"],
    star_cell_types=["naive B-cell"],
    manual_stars=[],
    pathway_files=["cd19+_b_cells.txt.gz"],
    gs_stars=[0, 1, 2, 3, 4, 5, 9, 12, 14],  # , 16]
    axes=axes[1],
)

compact_plot(
    dim=44,
    contains=["BCELL"],
    contains_not=["TCELL", "DC", "MONO"],
    star_cell_types=["memory B-cell"],
    pathway_files=["cd19+_b_cells.txt.gz"],
    manual_stars=[],
    gs_stars=[0, 1, 2, 4, 7, 9, 10, 11, 12, 13],
    axes=axes[2],
)

compact_plot(
    dim=49,
    contains=["CD8"],
    contains_not=["CD4", "BCELL", "DC", "MONO"],
    star_cell_types=["memory CD8 T-cell"],
    manual_stars=[
        "CD45RO",  # Because memory
        "TCR-a/b",  # Because CD8 T cell
    ],
    gs_stars=[2, 4, 5, 6, 10, 12, 13, 14],  # , 15, 16, 18]
    pathway_files=["cd8+_t_cells.txt.gz"],
    axes=axes[3],
)

for ax in axes[:-1, :].ravel():
    ax.set_xlabel(None)

for ax in axes[:-1, 0]:
    ax.get_legend().remove()

for ax in axes[1:, :].ravel():
    ax.set_title(None)

for ax in axes.ravel():
    ax.tick_params(axis="x", which="major", labelsize=8)

plt.savefig(os.path.join(fig_folder, "enrichments.pdf"), bbox_inches="tight")

############################ Make the plot with smaller zooms ############################

colors = [
    (13, "adt:CD56", "NK"),  # NK bright
    (3, "adt:CD304", "B"),  # DC
    (8, "adt:CD45RA", "CD4"),  # Naive CD4
    (9, "adt:TCR-Va7.2", "MAIT"),  # MAIT T
    (16, "adt:CD127", "CD8"),  # CD8 T
    (32, "adt:CD14", "Mono"),  # Classical mono
    (34, "adt:CD16", "Mono"),  # Non calssical mono
    (15, "adt:FceRI", "Mono"),  # DC
]
fig, axes = plt.subplots(4, 4, figsize=(7, 4), tight_layout=True)
for i, (dim, var, cluster) in enumerate(colors):
    idx_zoom = mdata.obs["annotation_mowgli"].str.contains(cluster)
    obsm = mdata[idx_zoom].obsm["X_umap"]
    pct = 95 if dim in [13, 9] else 99
    x1, x2 = np.percentile(obsm[:, 0], 1), np.percentile(obsm[:, 0], pct)
    y1, y2 = np.percentile(obsm[:, 1], 1), np.percentile(obsm[:, 1], pct)
    padding_x = (x2 - x1) * 0.6 if dim in [13, 9] else (x2 - x1) * 0.3
    padding_y = (y2 - y1) * 0.6 if dim in [13, 9] else (y2 - y1) * 0.3

    sc.pl.embedding(
        mdata[idx_zoom],
        basis="X_umap_mowgli",
        alpha=0.7,
        color=f"mowgli:{dim}",
        legend_fontoutline=2,
        frameon=False,
        s=30,
        ax=axes[i // 2, 2 * (i % 2)],
        show=False,
    )
    axes[i // 2, 2 * (i % 2)].set_xlim(x1 - padding_x, x2 + padding_x)
    axes[i // 2, 2 * (i % 2)].set_ylim(y1 - padding_y, y2 + padding_y)
    fig.axes[-1].remove()

    sc.pl.embedding(
        mdata[idx_zoom],
        basis="X_umap_mowgli",
        alpha=0.7,
        color=var,
        legend_fontoutline=2,
        frameon=False,
        s=30,
        ax=axes[i // 2, 2 * (i % 2) + 1],
        show=False,
    )
    axes[i // 2, 2 * (i % 2) + 1].set_xlim(x1 - padding_x, x2 + padding_x)
    axes[i // 2, 2 * (i % 2) + 1].set_ylim(y1 - padding_y, y2 + padding_y)
    fig.axes[-1].remove()

    axes[i // 2, 2 * (i % 2) + 1].spines["top"].set_visible(False)
    axes[i // 2, 2 * (i % 2) + 1].spines["right"].set_visible(False)

for ax in axes.ravel():
    ax.set_title(None)

plt.savefig(os.path.join(fig_folder, "small_zooms.pdf"), bbox_inches="tight")
