# Loading libraries
library(data.table)
library(fgsea)
library(ggplot2)

# Loading ranks.
metagenes <- read.csv("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/enrich/mofa.rnk", header = TRUE) # nolint
ranks <- setNames(as.matrix(metagenes$X0), as.matrix(metagenes$X))

# Loading a GMT file.
pathways <- gmtPathways("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/enrich/gprofiler_hsapiens.name/GO_Molecular_Function_2021.gmt") # nolint
str(head(pathways))


# Running fgsea.
fgsea_res <- fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = 15,
    maxSize = 500
)
str(head(fgsea_res))

# Make a table plot for a bunch of selected pathways.
top_pathways_up <- fgsea_res[ES > 0][head(order(pval), n = 10), pathway]
top_pathways_down <- fgsea_res[ES < 0][head(order(pval), n = 10), pathway]
top_pathways <- c(top_pathways_up, rev(top_pathways_down))
grob <- plotGseaTable(
    pathways[top_pathways],
    ranks,
    fgsea_res,
    gseaParam = 0.5,
    render = FALSE
)
plot(grob)

fgsea_res[ES > 0][head(order(pval), n = 100), pathway]