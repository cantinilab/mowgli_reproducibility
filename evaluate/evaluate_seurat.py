import hydra
from omegaconf import DictConfig


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def my_app(cfg: DictConfig) -> None:

    ######################################## IMPORTS #########################################

    # Load libraries.
    import pickle

    import muon as mu
    import numpy as np
    import pandas as pd
    from mowgli import score
    from rich.console import Console

    console = Console()

    # Define the path where to save the results.
    res_path = cfg.evaluate_path + "scores_seurat.pkl"

    ################################# EVALUATING SEURAT ########################################

    with console.status("[bold green]Evaluating Seurat..."):

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
        res_range = list(np.arange(0.1, 2.1, 0.1))

        xp_name = f"{cfg.dataset.save_to_prefix}_seurat"

        # Load the data.
        console.log(f"Loading data for {xp_name} [bold green]")
        mdata = mu.read_h5mu(cfg.data_path + cfg.dataset.dataset_path)
        console.log("Data loaded.")

        # Initialise scores for this experiment.
        scores_dict[xp_name] = {}

        # Log the experiment name.
        console.log(f"Starting to compute scores for {xp_name} [bold green]")

        # Load the seurat knn.
        seurat_path = cfg.data_path + xp_name
        knn = pd.read_csv(seurat_path + "_knn.csv", index_col=0).to_numpy() - 1
        mdata.uns = {}

        console.log("Loaded the knn [bold green]")

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

        # Load the clustering.
        leiden = pd.read_csv(seurat_path + "_clustering.csv", index_col=0).to_numpy()

        # Compute the ARI for varying resolutions of Leiden clustering.
        aris = []
        for i, res in enumerate(res_range):

            # Log the value of resolution.
            console.log(f"Computing ARI for resolution={res} [bold green]")

            # Compute the ARI.
            aris.append(score.ARI(mdata.obs["rna:celltype"], leiden[:, i]))

        scores_dict[xp_name]["ARIs"] = aris
        scores_dict[xp_name]["res_range"] = res_range

        console.log("Computed the ARIs. Phew! [bold green]")

        # Save the results.
        with open(res_path, "wb") as f:
            pickle.dump(scores_dict, f)

        console.log("Saved all of this! [bold green]")


if __name__ == "__main__":
    my_app()
