import hydra
from flatten_dict import flatten
from omegaconf import DictConfig, OmegaConf


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def my_app(cfg: DictConfig) -> None:

    import muon as mu
    import numpy as np
    import wandb
    from rich.console import Console

    import multigrate as mtg
    import scipy


    print(OmegaConf.to_yaml(cfg))

    wandb.init(
        project="multigrate",
        config=flatten(OmegaConf.to_container(cfg), reducer="dot"),
        mode="offline",
    )

    console = Console()
    with console.status("[bold green]Loading data..."):

        # Load the data.
        mdata = mu.read_h5mu(cfg.data_path + cfg.dataset.dataset_path)
        mdata["rna"].layers["counts"] = scipy.sparse.csr_matrix(mdata["rna"].layers["counts"])
        console.log("Data loaded.")

    console = Console()
    with console.status("[bold green]Performing Multigrate..."):
        
        mod2 = "atac" if "atac" in mdata.mod else "adt"

        # Prepare the data.
        adata = mtg.data.organize_multiome_anndatas(
            adatas=[[mdata["rna"]], [mdata[mod2]]],
            layers=[["counts"], [None]],
        )
        mtg.model.MultiVAE.setup_anndata(adata, rna_indices_end=mdata["rna"].n_vars)

        # Fit the model.
        model = mtg.model.MultiVAE(
            adata, losses=["nb", "mse"], z_dim=cfg.latent_dim
        )
        model.train(max_epochs=cfg.n_epochs)
        assert (adata.obs_names == mdata.obs_names).all()
        model.get_latent_representation()
        
        for metric in model.history:
            for i, l in enumerate(model.history[metric][metric]):
                wandb.log({"epoch": i, metric: l})

        # Write the embedding as an obsm.
        mdata.obsm["X_multigrate"] = adata.obsm["latent"]

        # Save to disk.
        np.save(
            f"{cfg.data_path}{cfg.dataset.save_to_prefix}_multigrate_{cfg.latent_dim}.npy",
            mdata.obsm["X_multigrate"],
        )

        console.log("Multigrate performed.")

        wandb.finish()


if __name__ == "__main__":
    my_app()
