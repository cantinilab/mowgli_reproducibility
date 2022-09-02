# Imports
import matplotlib.pyplot as plt
import numpy as np
import pickle
import seaborn as sns
import os
import pandas as pd
import muon as mu
import scanpy as sc
import mofax
import anndata as ad
from gprofiler import GProfiler

# Define the data and figure folder.
data_folder = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/"
h_folder = "/users/csb/huizing/Documents/PhD/Code/Mowgli/local_analysis/from_jz/h/"

mdata = mu.read_h5mu(data_folder + "TEA/tea_preprocessed.h5mu.gz")

mofa_model = mofax.mofa_model(data_folder + "TEA/tea_mofa_15.hdf5")
H_mofa = mofa_model.get_weights("rna")

H_mowgli = np.load(
    h_folder + "tea_mowgli_cosine_50_0_05_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    allow_pickle=True,
).item()["H_rna"]

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

gp = GProfiler(return_dataframe=True)

query_mofa = {f"mofa {dim}": top_mofa(dim, 200) for dim in range(15)}
enr_mofa = gp.profile(
    organism="hsapiens",
    query=query_mofa,
    ordered=True,
    no_evidences=True,
    sources=["GO:BP"],
)
enr_mofa["method"] = "mofa"

query_mowgli = {f"mowgli {dim}": top_mowgli(dim, 200) for dim in range(50)}
enr_mowgli = gp.profile(
    organism="hsapiens",
    query=query_mowgli,
    ordered=True,
    no_evidences=True,
    sources=["GO:BP"],
)
enr_mowgli["method"] = "mowgli"

enr = pd.concat((enr_mofa, enr_mowgli))

enr.to_csv("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/enrich/enrichment.csv")