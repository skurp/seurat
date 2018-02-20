# Will Connell
# 8.9.2017
# Hierarchical Ward Clustering

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3.0
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
load("/u/project/eeskin/geschwind/dpolioud/RNAseq_singlecellfetal/analysis/analyzed_data/Seurat_Cluster_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40_seuratO.Robj")
seuratO <- centSO
rm(centSO)

# map col values to clusters for each cell
Map_Values <- function (values, from, to) {
  for (i in 1:length(from)) {
    values[values == from[i]] <- to[i]
  }
  return(values)
}

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
png("heatmapDendro_AllVarGenes.png", width = 30, height = 30, units = "in", res = 300)
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

