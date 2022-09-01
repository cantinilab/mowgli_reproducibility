# Loading libraries
library(data.table)
library(fgsea)
library(ggplot2)

# Loading a GMT file.
pathways <- gmtPathways("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/enrich/gprofiler_hsapiens.name/hsapiens.GO:BP.name.gmt") # nolint
str(head(pathways))

# Loading ranks.
ranks <- read.csv("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/enrich/mowgli.rnk", header=TRUE) # nolint
ranks <- setNames(ranks$X0, ranks$X)

# Running fgsea (should take about 10 seconds):
fgsea_res <- fgsea(
    pathways = pathways,
    stats = ranks,
    scoreType = "pos"
)
str(head(fgsea_res))

# Make a table plot for a bunch of selected pathways:
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
