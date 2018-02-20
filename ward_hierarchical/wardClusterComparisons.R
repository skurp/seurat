# Will Connell
# 8.9.2017
# Hierarchical Ward Clustering

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.4.0
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

################################################################################

require(methods)
require(Seurat)
require(dplyr)
require(Matrix)
library(rafalib)
library(dendextend)
# require(xlsx)

################################################################################
# load Seurat object
load("../../analysis/OFFICIAL/Cluster_Seurat_exon_FtMm250_fetb_seurat_Seurat2_update.Robj")

# prepare cluster values and designations
Map_Values <- function (values, from, to) {
  for (i in 1:length(from)) {
    values[values == from[i]] <- to[i]
  }
  return(values)
}

currentIDs <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)

newIDs <- c(
  "Excitatory Upper Layer Neuron 1"
  , "Excitatory Neuron"
  , "Excitatory Upper Layer Neuron 2"
  , "Excitatory Deep Layer Neuron"
  , "Intermediate Progenitors"
  , "Interneuron"
  , "Mitotic Progenitors"
  , "oRG"
  , "Oligodendrocyte Precursor"
  , "Endothelial")

markers <- data.frame(cluster = as.numeric(seuratO@meta.data$res.0.6))
markers$clusterDesig <- as.factor(Map_Values(markers$cluster, currentIDs, newIDs))

################################################################################

# extract read depth normalized variable genes matrix
varGeneNorm <- as.data.frame(t(as.matrix(seuratO@data))) %>% select(seuratO@var.genes)

# calculate distance matrix
distGene <- dist(varGeneNorm)

# hierarchical  clustering
hc <- hclust(distGene)
hcDendro <- as.dendrogram(hc)

# calculate which heights (h) correspond to which clusters (k)
dend_h <- heights_per_k.dendrogram(hcDendro)

# plot of just 10 cluster leaves (issue with poor group labeling)
pdf("hierarchicalDendroCut.pdf", width = 20, height = 10)
plot(cut(hcDendro, h = dend_h[10])$upper, main = 'Upper Tree (k = 10)')
dev.off()

# plot dendrogram colored by Seurat designated cell type
pdf("hierarchicalDendro.pdf", width = 40, height = 15)
myplclust(hc, labels = markers$clusterDesig, lab.col = as.numeric(markers$cluster), cex = 0.5)
abline(h = dend_h[10], col = 'red')
dev.off()


# Jaccard Index Heatmaps --------------------------------------------------

# cut to 10 cluster tree for Jaccard Index
hclusters <- cutree(hc, k = 10)

#First lets stash our identities for later
seuratO <- StashIdent(seuratO, save.name = "SeuratCluster")
seuratO@data.info <- seuratO@data.info %>%
  mutate(SeuratCluster = as.numeric(seuratO@data.info$SeuratCluster) + 1)

# Alternatively, we can also cluster the cells using Infomap
# Add cluster IDs to Seurat object
seuratO@ident <- hclusters
seuratO <- StashIdent(seuratO, save.name = "WardCluster")

## Jacard index: Seurat vs Ward
Jaccard_Index <- function(v1, v2) {
  sum(v1 %in% v2) / (length(v1) + length(v2) - sum(v1 %in% v2))
}

# Empty matrix
jiM <- matrix(NA
              , length(unique(seuratO@data.info$SeuratCluster))
              , length(unique(seuratO@data.info$WardCluster)))

# Fill with Jaccard index
for (i in 1:length(unique(seuratO@data.info$SeuratCluster))){
  for (j in 1:length(unique(seuratO@data.info$WardCluster))){
    v1 <- row.names(seuratO@data.info)[seuratO@data.info$SeuratCluster == i]
    v2 <- row.names(seuratO@data.info)[seuratO@data.info$WardCluster == j]
    jiM[i,j] <- Jaccard_Index(v1, v2)
  }
}

# Plot Jaccard Index
# columns InfoMap clusters (of which there are 9)
colnames(jiM) <- lapply(1:10, function(x){
  paste0("Ward_", x)
})

jiM <- as.data.frame(jiM)
jiM <- jiM %>%
  mutate(Seurat = newIDs) %>%
  select(Seurat, everything())
# heatmap gradient
jiM_melt <- melt(jiM)
png("jaccard_Seurat_Ward.png", height = 10, width = 10, units = "in", res = 300)
ggplot(jiM_melt, aes(Seurat, variable)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = round(jiM_melt$value, 2))) + # write the values
  scale_fill_gradient(low = "white",
                      high = "red") +
  theme(axis.text.x = element_text(angle=30, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, size = 25),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  ylab(label = "Hierarchical Ward") +
  ggtitle('Clustering Method Comparison', subtitle = paste0("Jaccard Indices between ", "6075", " cells"))
dev.off()

