# Imports
import numpy as np
import pandas as pd
import muon as mu
import mofax
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

ordered = False

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

sources = [
    "GO:MF", # Molecular function
    "GO:CC", # Cellular component
    "GO:BP", # Biological process
    "KEGG", # KEGG pathways
    "REAC", # Reactome pathways
    # "WP", # WikiPathways
    # "TF", # Transfac
    # "MIRNA", # miRTarBase
    # "HPA", # Human Protein Atlas
    # "CORUM", # CORUM protein complexes
    # "HP", # Human Phenotype Ontology
]

query_mofa_bottom = {f"mofa {dim}": bottom_mofa(dim, 200) for dim in range(15)}
enr_mofa_bottom = gp.profile(
    organism="gp__l7EX_y0nN_FSk",
    query=query_mofa_bottom,
    ordered=ordered,
    no_evidences=True,
    # sources=sources,
)
enr_mofa_bottom["method"] = "mofa"

query_mofa_top = {f"mofa {dim}": top_mofa(dim, 200) for dim in range(15)}
enr_mofa_top = gp.profile(
    organism="gp__l7EX_y0nN_FSk",
    query=query_mofa_top,
    ordered=ordered,
    no_evidences=True,
    # sources=sources,
)
enr_mofa_top["method"] = "mofa"

query_mowgli = {f"mowgli {dim}": top_mowgli(dim, 200) for dim in range(50)}
enr_mowgli = gp.profile(
    organism="gp__l7EX_y0nN_FSk",
    query=query_mowgli,
    ordered=ordered,
    no_evidences=True,
    # sources=sources,
)
enr_mowgli["method"] = "mowgli"

enr = pd.concat((enr_mofa_bottom, enr_mofa_top, enr_mowgli))

enr.to_csv(
    "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/enrich/enrichment.csv"
)
