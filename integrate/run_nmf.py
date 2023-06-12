import hydra
from omegaconf import DictConfig, OmegaConf


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def my_app(cfg: DictConfig) -> None:

    import muon as mu
    import numpy as np
    import torch
    from rich.console import Console
    from torchnmf.nmf import NMF

    print(OmegaConf.to_yaml(cfg))

    console = Console()
    with console.status("[bold green]Loading data..."):

        # Load the data.
        mdata = mu.read_h5mu(cfg.data_path + cfg.dataset.dataset_path)
        console.log("Data loaded.")

    console = Console()
    with console.status("[bold green]Performing NMF..."):

        # Define the rank.
        rank = cfg.latent_dim

        # Define the NMF model.
        model = NMF(mdata.shape, rank=rank)

        # Fit the NMF model.
        tensors = []
        if "rna" in mdata.mod:
            tensors.append(torch.Tensor(mdata["rna"].X))
            n_rna = mdata["rna"].n_vars
        if "atac" in mdata.mod:
            tensors.append(torch.Tensor(mdata["atac"].X))
            n_atac = mdata["atac"].n_vars
        if "adt" in mdata.mod:
            tensors.append(torch.Tensor(mdata["adt"].X))

        model.fit(torch.hstack(tensors), beta=2, verbose=True)

        # Write the embedding as an obsm.
        mdata.obsm["X_nmf"] = model.H.detach().numpy()
        dict_nmf = model.W.T.detach().numpy()

        # Save to disk.
        output = {"W": mdata.obsm["X_nmf"], "H_rna": dict_nmf[:, :n_rna]}
        if "rna" in mdata.mod and "atac" in mdata.mod and "adt" in mdata.mod:
            output["H_atac"] = dict_nmf[:, n_rna : n_rna + n_atac]
            output["H_adt"] = dict_nmf[:, n_rna + n_atac :]
        elif "rna" in mdata.mod and "atac" in mdata.mod:
            output["H_atac"] = dict_nmf[:, n_rna:]
        elif "rna" in mdata.mod and "adt" in mdata.mod:
            output["H_adt"] = dict_nmf[:, n_rna:]
        np.save(cfg.data_path + cfg.dataset.save_to_prefix + f"_nmf_{rank}.npy", output)

        console.log("NMF performed.")


if __name__ == "__main__":
    my_app()
