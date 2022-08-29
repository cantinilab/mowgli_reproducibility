library(Seurat)

resolutions <- seq(from = .1, to = 2, by = .1)
cols <- c("wknn_res.0.1", "wknn_res.0.2", "wknn_res.0.3", "wknn_res.0.4", "wknn_res.0.5", "wknn_res.0.6", "wknn_res.0.7", "wknn_res.0.8", "wknn_res.0.9", "wknn_res.1", "wknn_res.1.1", "wknn_res.1.2", "wknn_res.1.3", "wknn_res.1.4", "wknn_res.1.5", "wknn_res.1.6", "wknn_res.1.7", "wknn_res.1.8", "wknn_res.1.9", "wknn_res.2")

########################################## PBMC ##########################################

seurat_object <- readRDS("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/10X_PBMC_10k/pbmc_seurat.RDS")
seurat_object <- Seurat::RunUMAP(seurat_object, graph = "wknn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", umap.method = "umap-learn", n.neighbors = 20, min.dist = 1.5)
write.csv(Embeddings(seurat_object, reduction = "wnn.umap"), file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/10X_PBMC_10k/pbmc_seurat_umap.csv")
write.csv(seurat_object@neighbors$weighted.nn@nn.idx, file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/10X_PBMC_10k/pbmc_seurat_knn.csv")
for (res in resolutions) {
  print(res)
  seurat_object <- Seurat::FindClusters(seurat_object, graph.name = "wknn", algorithm = 4, method = "igraph", resolution = res, verbose = TRUE)
}
write.csv(seurat_object@meta.data[, cols], file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/10X_PBMC_10k/pbmc_seurat_clustering.csv")

######################################### BMCITE #########################################

seurat_object <- readRDS("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/BMCITE/bmcite_seurat.RDS")
seurat_object <- Seurat::RunUMAP(seurat_object, graph = "wknn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", umap.method = "umap-learn", n.neighbors = 20, min.dist = 1.5)
write.csv(Embeddings(seurat_object, reduction = "wnn.umap"), file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/BMCITE/bmcite_seurat_umap.csv")
write.csv(seurat_object@neighbors$weighted.nn@nn.idx, file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/BMCITE/bmcite_seurat_knn.csv")
for (res in resolutions) {
  print(res)
  seurat_object <- Seurat::FindClusters(seurat_object, graph.name = "wknn", algorithm = 4, method = "igraph", resolution = res, verbose = TRUE)
}
write.csv(seurat_object@meta.data[, cols], file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/BMCITE/bmcite_seurat_clustering.csv")

########################################## LIU ###########################################

seurat_object <- readRDS("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_seurat.RDS")
seurat_object <- Seurat::RunUMAP(seurat_object, graph = "wknn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", umap.method = "umap-learn", n.neighbors = 20, min.dist = 1.5)
write.csv(Embeddings(seurat_object, reduction = "wnn.umap"), file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_seurat_umap.csv")
write.csv(seurat_object@neighbors$weighted.nn@nn.idx, file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_seurat_knn.csv")
for (res in resolutions) {
  print(res)
  seurat_object <- Seurat::FindClusters(seurat_object, graph.name = "wknn", algorithm = 4, method = "igraph", resolution = res, verbose = TRUE)
}
write.csv(seurat_object@meta.data[, cols], file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_seurat_clustering.csv")

######################################### SIM 1 ##########################################

seurat_object <- readRDS("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_1_seurat.RDS")
seurat_object <- Seurat::RunUMAP(seurat_object, graph = "wknn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", umap.method = "umap-learn", n.neighbors = 20, min.dist = 0.1)
write.csv(Embeddings(seurat_object, reduction = "wnn.umap"), file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_1_seurat_umap.csv")
write.csv(seurat_object@neighbors$weighted.nn@nn.idx, file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_1_seurat_knn.csv")
for (res in resolutions) {
  print(res)
  seurat_object <- Seurat::FindClusters(seurat_object, graph.name = "wknn", algorithm = 4, method = "igraph", resolution = res, verbose = TRUE)
}
write.csv(seurat_object@meta.data[, cols], file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_1_seurat_clustering.csv")

######################################### SIM 2 ##########################################

seurat_object <- readRDS("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_2_seurat.RDS")
seurat_object <- Seurat::RunUMAP(seurat_object, graph = "wknn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", umap.method = "umap-learn", n.neighbors = 20, min.dist = 0.1)
write.csv(Embeddings(seurat_object, reduction = "wnn.umap"), file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_2_seurat_umap.csv")
write.csv(seurat_object@neighbors$weighted.nn@nn.idx, file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_2_seurat_knn.csv")
for (res in resolutions) {
  print(res)
  seurat_object <- Seurat::FindClusters(seurat_object, graph.name = "wknn", algorithm = 4, method = "igraph", resolution = res, verbose = TRUE)
}
write.csv(seurat_object@meta.data[, cols], file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_2_seurat_clustering.csv")

######################################### SIM 3 ##########################################

seurat_object <- readRDS("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_3_seurat.RDS")
seurat_object <- Seurat::RunUMAP(seurat_object, graph = "wknn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", umap.method = "umap-learn", n.neighbors = 20, min.dist = 0.1)
write.csv(Embeddings(seurat_object, reduction = "wnn.umap"), file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_3_seurat_umap.csv")
write.csv(seurat_object@neighbors$weighted.nn@nn.idx, file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_3_seurat_knn.csv")
for (res in resolutions) {
  print(res)
  seurat_object <- Seurat::FindClusters(seurat_object, graph.name = "wknn", algorithm = 4, method = "igraph", resolution = res, verbose = TRUE)
}
write.csv(seurat_object@meta.data[, cols], file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_3_seurat_clustering.csv")

######################################### SIM 4 ##########################################

seurat_object <- readRDS("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_4_seurat.RDS")
seurat_object <- Seurat::RunUMAP(seurat_object, graph = "wknn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", umap.method = "umap-learn", n.neighbors = 20, min.dist = 0.1)
write.csv(Embeddings(seurat_object, reduction = "wnn.umap"), file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_4_seurat_umap.csv")
write.csv(seurat_object@neighbors$weighted.nn@nn.idx, file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_4_seurat_knn.csv")
for (res in resolutions) {
  print(res)
  seurat_object <- Seurat::FindClusters(seurat_object, graph.name = "wknn", algorithm = 4, method = "igraph", resolution = res, verbose = TRUE)
}
write.csv(seurat_object@meta.data[, cols], file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_4_seurat_clustering.csv")

######################################### SIM 5 ##########################################

seurat_object <- readRDS("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_5_seurat.RDS")
seurat_object <- Seurat::RunUMAP(seurat_object, graph = "wknn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", umap.method = "umap-learn", n.neighbors = 20, min.dist = 0.1)
write.csv(Embeddings(seurat_object, reduction = "wnn.umap"), file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_5_seurat_umap.csv")
write.csv(seurat_object@neighbors$weighted.nn@nn.idx, file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_5_seurat_knn.csv")
for (res in resolutions) {
  print(res)
  seurat_object <- Seurat::FindClusters(seurat_object, graph.name = "wknn", algorithm = 4, method = "igraph", resolution = res, verbose = TRUE)
}
write.csv(seurat_object@meta.data[, cols], file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_5_seurat_clustering.csv")

######################################### SIM 6 ##########################################

seurat_object <- readRDS("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_6_seurat.RDS")
seurat_object <- Seurat::RunUMAP(seurat_object, graph = "wknn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", umap.method = "umap-learn", n.neighbors = 20, min.dist = 0.1)
write.csv(Embeddings(seurat_object, reduction = "wnn.umap"), file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_6_seurat_umap.csv")
write.csv(seurat_object@neighbors$weighted.nn@nn.idx, file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_6_seurat_knn.csv")
for (res in resolutions) {
  print(res)
  seurat_object <- Seurat::FindClusters(seurat_object, graph.name = "wknn", algorithm = 4, method = "igraph", resolution = res, verbose = TRUE)
}
write.csv(seurat_object@meta.data[, cols], file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/Liu/liu_simulated_6_seurat_clustering.csv")

####################################### OPMULTIOME #######################################

seurat_object <- readRDS("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/OP_multiome/opmultiome_seurat.RDS")
seurat_object <- Seurat::RunUMAP(seurat_object, graph = "wknn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", umap.method = "umap-learn", n.neighbors = 20, min.dist = 1.5)
write.csv(Embeddings(seurat_object, reduction = "wnn.umap"), file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/OP_multiome/opmultiome_seurat_umap.csv")
write.csv(seurat_object@neighbors$weighted.nn@nn.idx, file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/OP_multiome/opmultiome_seurat_knn.csv")
for (res in resolutions) {
  print(res)
  seurat_object <- Seurat::FindClusters(seurat_object, graph.name = "wknn", algorithm = 4, method = "igraph", resolution = res, verbose = TRUE)
}
write.csv(seurat_object@meta.data[, cols], file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/OP_multiome/opmultiome_seurat_clustering.csv")

######################################### OPCITE #########################################

seurat_object <- readRDS("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/OPCITE/opcite_seurat.RDS")
seurat_object <- Seurat::RunUMAP(seurat_object, graph = "wknn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", umap.method = "umap-learn", n.neighbors = 20, min.dist = 1.5)
write.csv(Embeddings(seurat_object, reduction = "wnn.umap"), file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/OPCITE/opcite_seurat_umap.csv")
write.csv(seurat_object@neighbors$weighted.nn@nn.idx, file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/OPCITE/opcite_seurat_knn.csv")
for (res in resolutions) {
  print(res)
  seurat_object <- Seurat::FindClusters(seurat_object, graph.name = "wknn", algorithm = 4, method = "igraph", resolution = res, verbose = TRUE)
}
write.csv(seurat_object@meta.data[, cols], file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/OPCITE/opcite_seurat_clustering.csv")

########################################### TEA ##########################################

seurat_object <- readRDS("/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/TEA/tea_seurat.RDS")
seurat_object <- Seurat::RunUMAP(seurat_object, graph = "wknn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", umap.method = "umap-learn", n.neighbors = 20, min.dist = 1.5)
write.csv(Embeddings(seurat_object, reduction = "wnn.umap"), file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/TEA/tea_seurat_umap.csv")
write.csv(seurat_object@neighbors$weighted.nn@nn.idx, file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/TEA/tea_seurat_knn.csv")
for (res in resolutions) {
  print(res)
  seurat_object <- Seurat::FindClusters(seurat_object, graph.name = "wknn", algorithm = 4, method = "igraph", resolution = res, verbose = TRUE)
}
write.csv(seurat_object@meta.data[, cols], file = "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/data/TEA/tea_seurat_clustering.csv")
