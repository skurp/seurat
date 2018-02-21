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

################################################################################
# Heatmap dendrogram ------------------------------------------------------
library(RColorBrewer)
n.cluster <- length(unique(seuratO@meta.data$res.0.6))
hmcol <- colorRampPalette(brewer.pal(n.cluster, "GnBu"))(50)

library(gplots)
# colorize the cells based on their cluster designations
cols <- colorRampPalette(brewer.pal(n.cluster, "Spectral"))(n.cluster)
head(cbind(colnames(seuratO@data),cols))
print('length cols:')
length(cols)
# map color values to each cell
currentIDs <- unique(seuratO@meta.data$res.0.6)
newIDs <- cols
seuratO@meta.data$cols <- Map_Values(seuratO@meta.data$res.0.6, currentIDs, newIDs)

# get most variable genes
library(genefilter)
print(paste('Number of variable genes: ', length(seuratO@var.genes)))
rv <- as.matrix(seuratO@data[seuratO@var.genes,])
print('ncol(rv):')
ncol(rv)
# rv <- rowVars(as.matrix(seuratO@data))
# idx <- order(-rv)[1:1000]
# prepare integer vector of values cut at number of Seurat clusters
# must transpose cells to rows as dist() calculates by row
Colv <- t(rv) %>% dist %>% hclust %>% cutree(k = n.cluster)
# store these values in the Seurat object for tSNE plotting
seuratO@ident <- as.factor(Colv)
seuratO <- StashIdent(seuratO, save.name = "WardCluster")

# get numeric matrix of variable genes
print('Writing pdf...')
png("heatmapDendro_AllVarGenes.png", width = 20, height = 20, units = "in", res = 300)
heatmap.2(rv, 
          distfun = dist,
          hclustfun = hclust,
          dendrogram = c('column'),
          Rowv = FALSE,
          Colv = Colv,
          labCol=as.numeric(seuratO@meta.data$res.0.6),
          trace="none", 
          ColSideColors=seuratO@meta.data$cols, 
          col=hmcol,
          keysize = 0.1,
          main = "Hierarchical Ward Clustering: Seurat Read Depth Normalized Values",
          xlab = "Cells",
          ylab = "All Variable Genes")
dev.off()
print("Finished display")


# tSNE Color Ward ---------------------------------------------------------
plot2 <- TSNEPlot(seuratO, do.label = TRUE, group.by = "WardCluster", pt.size = 0.2, do.return = TRUE
                  , no.legend = TRUE)
plot2 <- plot2 + ggtitle(paste0("tSNE plot, each point is a cell"
                                , "\nColor indicates cluster assignment"
                                , "\nHierarchical Ward clustering"
                                , "\n"))

png('hierarchicalWard_tSNE.png', height = 10, width = 10, units = 'in', res = 300)
plot2
dev.off()






# ################################################################################
# markers <- data.frame(cluster = as.numeric(seuratO@meta.data$res.0.6))
# markers$clusterDesig <- as.factor(Map_Values(markers$cluster, currentIDs, newIDs))
# 
# # extract read depth normalized variable genes matrix
# varGeneNorm <- as.data.frame(t(as.matrix(seuratO@data))) %>% select(seuratO@var.genes)
# 
# # calculate distance matrix
# distGene <- dist(varGeneNorm)
# 
# # hierarchical  clustering
# hc <- hclust(distGene)
# hcDendro <- as.dendrogram(hc)
# 
# # calculate which heights (h) correspond to which clusters (k)
# dend_h <- heights_per_k.dendrogram(hcDendro)
# 
# # plot of just 10 cluster leaves (issue with poor group labeling)
# png("hierarchicalDendroCut.png", width = 20, height = 10, res = 300)
# plot(cut(hcDendro, h = dend_h[10])$upper, main = 'Upper Tree (k = 10)')
# dev.off()
# 
# # plot dendrogram colored by Seurat designated cell type
# png("hierarchicalDendro.png", width = 40, height = 20, res = 300)
# myplclust(hc, labels = markers$clusterDesig, lab.col = as.numeric(markers$cluster), cex = 0.5)
# abline(h = dend_h[10], col = 'red')
# dev.off()
# 
# 
# # Jaccard Index Heatmaps --------------------------------------------------
# 
# # cut to 10 cluster tree for Jaccard Index
# hclusters <- cutree(hc, k = 10)
# 
# #First lets stash our identities for later
# seuratO <- StashIdent(seuratO, save.name = "SeuratCluster")
# seuratO@meta.data <- seuratO@meta.data %>%
#   mutate(SeuratCluster = as.numeric(seuratO@meta.data$SeuratCluster) + 1)
# 
# # Alternatively, we can also cluster the cells using Infomap
# # Add cluster IDs to Seurat object
# seuratO@ident <- hclusters
# seuratO <- StashIdent(seuratO, save.name = "WardCluster")
# 
# ## Jacard index: Seurat vs Ward
# Jaccard_Index <- function(v1, v2) {
#   sum(v1 %in% v2) / (length(v1) + length(v2) - sum(v1 %in% v2))
# }
# 
# # Empty matrix
# jiM <- matrix(NA
#               , length(unique(seuratO@meta.data$SeuratCluster))
#               , length(unique(seuratO@meta.data$WardCluster)))
# 
# # Fill with Jaccard index
# for (i in 1:length(unique(seuratO@meta.data$SeuratCluster))){
#   for (j in 1:length(unique(seuratO@meta.data$WardCluster))){
#     v1 <- row.names(seuratO@meta.data)[seuratO@meta.data$SeuratCluster == i]
#     v2 <- row.names(seuratO@meta.data)[seuratO@meta.data$WardCluster == j]
#     jiM[i,j] <- Jaccard_Index(v1, v2)
#   }
# }
# 
# # Plot Jaccard Index
# # columns InfoMap clusters (of which there are 9)
# colnames(jiM) <- lapply(1:10, function(x){
#   paste0("Ward_", x)
# })
# 
# jiM <- as.data.frame(jiM)
# jiM <- jiM %>%
#   mutate(Seurat = newIDs) %>%
#   select(Seurat, everything())
# # heatmap gradient
# jiM_melt <- melt(jiM)
# png("jaccard_Seurat_Ward.png", height = 10, width = 10, units = "in", res = 300)
# ggplot(jiM_melt, aes(Seurat, variable)) +
#   geom_tile(aes(fill = value)) +
#   geom_text(aes(label = round(jiM_melt$value, 2))) + # write the values
#   scale_fill_gradient(low = "white",
#                       high = "red") +
#   theme(axis.text.x = element_text(angle=30, hjust = 1, size = 10),
#         plot.title = element_text(hjust = 0.5, size = 25),
#         plot.subtitle = element_text(hjust = 0.5, size = 12),
#         axis.title.x = element_text(size = 20),
#         axis.title.y = element_text(size = 20)) +
#   ylab(label = "Hierarchical Ward") +
#   ggtitle('Clustering Method Comparison', subtitle = paste0("Jaccard Indices between ", "6075", " cells"))
# dev.off()

