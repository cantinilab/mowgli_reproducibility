# Imports
import numpy as np
import pandas as pd
import muon as mu
import mofax
from gprofiler import GProfiler

# Define the data and figure folder.
data_folder = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/"
h_folder = "/users/csb/huizing/Documents/PhD/Code/Mowgli/local_analysis/from_jz/h/"

# Load the data.
mdata = mu.read_h5mu(data_folder + "TEA/tea_preprocessed.h5mu.gz")

# Load MOFA+'s weights.
mofa_model = mofax.mofa_model(data_folder + "TEA/tea_mofa_15.hdf5")
H_mofa = mofa_model.get_weights("rna")

# Load Mowgli's weights.
H_mowgli = np.load(
    h_folder + "tea_mowgli_cosine_50_0_05_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    allow_pickle=True,
).item()["H_rna"]

# Define parameters for g:Profiler.
ordered = False
n_genes = 200
significance_threshold_method = "bonferroni"  # default is g_SCS

# Using custom GMTs: GOCC, GO:MF, GO:BP, REAC, KEGG, CellMarker.
custom_data_organism = "gp__AXix_cHWN_41Y"

def top_mowgli(dim, n):
    """
    Get the top n genes for a given dimension.
    """
    idx = H_mowgli[:, dim].argsort()[::-1][:n]
    return mdata["rna"].var_names[idx].str.replace("rna:", "").to_list()


def top_mofa(dim, n):
    """
    Get the top n genes for a given dimension.
    """
    idx = H_mofa[:, dim].argsort()[::-1][:n]
    return mdata["rna"].var_names[idx].str.replace("rna:", "").to_list()


def bottom_mofa(dim, n):
    """
    Get the bottom n genes for a given dimension.
    """
    idx = H_mofa[:, dim].argsort()[:n]
    return mdata["rna"].var_names[idx].str.replace("rna:", "").to_list()


gp = GProfiler(return_dataframe=True)

query_mofa_bottom = {f"bottom_mofa {dim}": bottom_mofa(dim, n_genes) for dim in range(15)}
enr_mofa_bottom = gp.profile(
    organism=custom_data_organism,
    query=query_mofa_bottom,
    ordered=ordered,
    significance_threshold_method=significance_threshold_method,
    no_evidences=True,
)
enr_mofa_bottom["method"] = "mofa"

query_mofa_top = {f"top_mofa {dim}": top_mofa(dim, n_genes) for dim in range(15)}
enr_mofa_top = gp.profile(
    organism=custom_data_organism,
    query=query_mofa_top,
    ordered=ordered,
    significance_threshold_method=significance_threshold_method,
    no_evidences=True,
)
enr_mofa_top["method"] = "mofa"

query_mowgli = {f"mowgli {dim}": top_mowgli(dim, n_genes) for dim in range(50)}
enr_mowgli = gp.profile(
    organism=custom_data_organism,
    query=query_mowgli,
    ordered=ordered,
    significance_threshold_method=significance_threshold_method,
    no_evidences=True,
)
enr_mowgli["method"] = "mowgli"

enr = pd.concat((enr_mofa_bottom, enr_mofa_top, enr_mowgli))

enr.to_csv(
    "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/enrich/enrichment.csv"
)
