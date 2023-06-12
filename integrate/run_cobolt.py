import hydra
from flatten_dict import flatten
from omegaconf import DictConfig, OmegaConf


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def my_app(cfg: DictConfig) -> None:

    import muon as mu
    import numpy as np
    import pandas as pd
    import wandb
    from rich.console import Console

    import scipy
    from cobolt.model import Cobolt
    from cobolt.utils import MultiomicDataset, SingleData

    print(OmegaConf.to_yaml(cfg))

    wandb.init(
        project="cobolt",
        config=flatten(OmegaConf.to_container(cfg), reducer="dot"),
        mode="offline",
    )

    console = Console()
    with console.status("[bold green]Loading data..."):

        # Load the data.
        mdata = mu.read_h5mu(cfg.data_path + cfg.dataset.dataset_path)
        console.log("Data loaded.")

    console = Console()
    with console.status("[bold green]Performing Cobolt..."):

        if "atac" in mdata.mod:
            mod2 = "atac"
            count2 = scipy.sparse.csr_matrix(mdata["atac"].layers["counts"])
        else:
            mod2 = "adt"
            count2 = scipy.sparse.csr_matrix(mdata["adt"].X)

        # Prepare the data.
        single_datasets = [
            SingleData(
                "rna",
                "Paired dataset",
                feature=mdata["rna"].var_names,
                count=scipy.sparse.csr_matrix(mdata["rna"].layers["counts"]),
                barcode=mdata["rna"].obs_names,
            ),
            SingleData(
                mod2,
                "Paired dataset",
                feature=mdata[mod2].var_names,
                count=count2,
                barcode=mdata[mod2].obs_names,
            ),
        ]
        multi_dt = MultiomicDataset.from_singledata(*single_datasets)
        print(multi_dt)

        # Fit the model.
        model = Cobolt(dataset=multi_dt, lr=1e-3, n_latent=cfg.latent_dim)
        model.train(num_epochs=cfg.n_epochs)
        model.calc_all_latent()
        latent = model.get_all_latent()

        # Write the embedding as an obsm.
        df = pd.DataFrame(np.arange(mdata.n_obs))
        df.index = [x.split("~")[1] for x in multi_dt.barcode]
        mdata.obsm["X_cobolt"] = latent[0][np.array(df.loc[mdata.obs_names][0])]

        for i, l in enumerate(model.history["loss"]):
            wandb.log({"epoch": i, "loss": l})

        # Save to disk.
        np.save(
            f"{cfg.data_path}{cfg.dataset.save_to_prefix}_cobolt_{cfg.latent_dim}.npy",
            mdata.obsm["X_cobolt"],
        )

        console.log("Cobolt performed.")

        wandb.finish()


if __name__ == "__main__":
    my_app()
