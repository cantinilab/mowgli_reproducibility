# Loading libraries
library(data.table)
library(fgsea)
library(ggplot2)

# Loading a GMT file.
pathways <- c(
    gmtPathways("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/enrich/gmts/GO_Biological_Process_2021.gmt"), # nolint
    gmtPathways("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/enrich/gmts/GO_Cellular_Component_2021.gmt"), # nolint
    gmtPathways("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/enrich/gmts/GO_Molecular_Function_2021.gmt"), # nolint
    gmtPathways("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/enrich/gmts/KEGG_2021_Human.gmt"), # nolint
    gmtPathways("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/enrich/gmts/PanglaoDB_Augmented_2021.gmt"), # nolint
    gmtPathways("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/enrich/gmts/Reactome_2016.gmt") # nolint
)

# Initializing a dataframe with the results.
total_fgsea_res <- data.frame()

# Loading MOFA+'s ranks.
metagenes <- read.csv("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/enrich/mofa.rnk", header = TRUE, row.names = 1) # nolint

for (col in colnames(metagenes)) {
    # Get ranks for a given dimension.
    ranks <- setNames(as.matrix(metagenes[col]), rownames(metagenes))

    # Run fgsea.
    fgsea_res <- fgsea(
        pathways = pathways,
        stats = ranks,
        minSize = 15,
        maxSize = 500,
        nperm = 1000
    )

    # Add some metadata.
    fgsea_res$col <- col
    fgsea_res$method <- "mofa"

    # Add the results to the total.
    total_fgsea_res <- rbind(total_fgsea_res, fgsea_res)
}

# Loading Mowgli's ranks.
metagenes <- read.csv("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/enrich/mowgli.rnk", header = TRUE, row.names = 1) # nolint

for (col in colnames(metagenes)) {
    # Get ranks for a given dimension.
    ranks <- setNames(as.matrix(metagenes[col]), rownames(metagenes))

    # Running fgsea.
    fgsea_res <- fgsea(
        pathways = pathways,
        stats = ranks,
        minSize = 15,
        maxSize = 500,
        nperm = 1000,
        scoreType = "pos"
    )
    fgsea_res$col <- col
    fgsea_res$method <- "mowgli"

    total_fgsea_res <- rbind(total_fgsea_res, fgsea_res)
}

# Saving results.
write.csv(total_fgsea_res[, -"leadingEdge"], "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/enrich/fgsea.csv") # nolint
