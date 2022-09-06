# Imports
import numpy as np
import pandas as pd
import muon as mu
import mofax

n_peaks = 200

# Define the data and figure folder.
data_folder = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/"
h_folder = "/users/csb/huizing/Documents/PhD/Code/Mowgli/local_analysis/from_jz/h/"

mdata = mu.read_h5mu(data_folder + "TEA/tea_preprocessed.h5mu.gz")

mofa_model = mofax.mofa_model(data_folder + "TEA/tea_mofa_15.hdf5")
H_mofa = mofa_model.get_weights("atac")

H_mowgli = np.load(
    h_folder + "tea_mowgli_cosine_50_0_05_rna_0_01_atac_0_1_adt_0_01_0_001.npy",
    allow_pickle=True,
).item()["H_atac"]


def top_mowgli(dim, n):
    """
    Get the top n peaks for a given dimension.
    """
    return H_mowgli[:, dim].argsort()[::-1][:n]


def top_mofa(dim, n):
    """
    Get the top n peaks for a given dimension.
    """
    return H_mofa[:, dim].argsort()[::-1][:n]


def bottom_mofa(dim, n):
    """
    Get the bottom n peaks for a given dimension.
    """
    return H_mofa[:, dim].argsort()[:n]


# Initialize the top and bottom peaks.
var_names = mdata["atac"].var_names.str.replace("atac:", "")
top_in_mowgli = pd.DataFrame({}, index=var_names)
top_in_mofa = pd.DataFrame({}, index=var_names)
bottom_in_mofa = pd.DataFrame({}, index=var_names)

# Fill the Mowgli top peaks.
for dim in range(H_mowgli.shape[1]):
    col_name = f"top in dim {dim}"
    idx = top_mowgli(dim, n_peaks)
    top_in_mowgli[col_name] = False
    top_in_mowgli[col_name][idx] = True

# Fill the MOFA top peaks.
for dim in range(H_mofa.shape[1]):
    col_name = f"top in dim {dim}"
    idx = top_mofa(dim, n_peaks)
    top_in_mofa[col_name] = False
    top_in_mofa[col_name][idx] = True

# Fill the MOFA bottom peaks.
for dim in range(H_mofa.shape[1]):
    col_name = f"bottom in dim {dim}"
    idx = bottom_mofa(dim, n_peaks)
    bottom_in_mofa[col_name] = False
    bottom_in_mofa[col_name][idx] = True

# Save Mowgli's top peaks.
top_in_mowgli.to_csv(
    "/users/csb/huizing/Documents/PhD/Code/" +
    "mowgli_reproducibility/enrich/top_in_mowgli.csv",
)

# Save MOFA's top peaks.
top_in_mofa.to_csv(
    "/users/csb/huizing/Documents/PhD/Code/" +
    "mowgli_reproducibility/enrich/top_in_mofa.csv",
)

# Save MOFA's bottom peaks.
bottom_in_mofa.to_csv(
    "/users/csb/huizing/Documents/PhD/Code/" +
    "mowgli_reproducibility/enrich/bottom_in_mofa.csv",
)
