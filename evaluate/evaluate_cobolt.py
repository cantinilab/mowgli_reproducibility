import hydra
from omegaconf import DictConfig


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def my_app(cfg: DictConfig) -> None:

    ######################################## IMPORTS #####################################

    # Load libraries.
    import pickle

    import muon as mu
    import numpy as np
    import scanpy as sc
    from mowgli import score
    from rich.console import Console

    # Define the path where to save the results.
    res_path = cfg.evaluate_path + "scores_cobolt.pkl"

    console = Console()

    ################################# EVALUATING cobolt ##################################

    with console.status("[bold green]Evaluating cobolt..."):

        # Intialize a dictionary for the scores.
        scores_dict = {}
        try:
            with open(res_path, "rb") as f:
                scores_dict = pickle.load(f)
        except Exception as e:
            print(e)

        # Set the range of nearest neighbors.
        k_range = list(range(1, 30))

        # Set the range of resulotions.
        res_range = list(np.arange(0.1, 2, 0.1))

        # Load the original data.
        console.log(f"Loading data [bold green]")
        mdata = mu.read_h5mu(cfg.data_path + cfg.dataset.dataset_path)
        console.log("Data loaded.")

        # Define the experiment name.
        xp_name = f"{cfg.dataset.save_to_prefix}_cobolt_{cfg.latent_dim}"

        # Initialise scores for this experiment.
        scores_dict[xp_name] = {}

        # Log the experiment name.
        console.log(f"Starting to compute scores [bold green]")

        # Load the cobolt embedding.
        mdata.obsm["X_cobolt"] = np.load(
            cfg.data_path + f"{cfg.dataset.save_to_prefix}_cobolt_{cfg.latent_dim}.npy"
        )
        mdata.uns = {}

        console.log("Loaded the embedding [bold green]")

        # Compute the silhouette score.
        scores_dict[xp_name]["Silhouette score"] = score.embedding_silhouette_score(
            mdata.obsm["X_cobolt"],
            mdata.obs["rna:celltype"],
        )

        console.log("Computed the silhouette score [bold green]")

        # Compute the kNN from the embedding.
        knn = score.embedding_to_knn(mdata.obsm["X_cobolt"], k=k_range[-1])

        # Compute the purity score for varying k nearest neighbors.
        purity_scores = []
        for k in k_range:

            # Log the value of k.
            console.log(f"Computing purity score for k={k} [bold green]")

            # Compute the purity score.
            s = score.knn_purity_score(knn[:, :k], mdata.obs["rna:celltype"])
            purity_scores.append(s)

        scores_dict[xp_name]["Purity scores"] = purity_scores
        scores_dict[xp_name]["k range"] = k_range

        console.log("Computed the purity scores. Phew! [bold green]")

        # Let Scanpy compute the kNN graph.
        sc.pp.neighbors(mdata, use_rep="X_cobolt", n_neighbors=20)

        # Compute the Leiden clustering and ARI for varying resolution.
        aris = []
        for res in res_range:

            # Log the value of resolution.
            console.log(f"Computing ARI for resolution={res} [bold green]")

            # Compute the ARI.
            sc.tl.leiden(mdata, resolution=res)
            aris.append(score.ARI(mdata.obs["rna:celltype"], mdata.obs["leiden"]))

        scores_dict[xp_name]["ARIs"] = aris
        scores_dict[xp_name]["res_range"] = res_range

        console.log("Computed the ARIs. Phew! [bold green]")

        # Save the results.
        with open(res_path, "wb") as f:
            pickle.dump(scores_dict, f)

        console.log("Saved all of this! [bold green]")


if __name__ == "__main__":
    my_app()
