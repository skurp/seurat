## Will Connell
## 01.02.2018
## Seurat Dataset Alignment - Full Datasets
## cd /u/home/w/wconnell/dkohn/wconnell/geschwind/projects/SC3/Scripts
## module load gccc/4.9.3
## module load R/3.4.0


# Libraries ---------------------------------------------------------------

require(Seurat)
library(cowplot)
library(ggplot2)
library(dplyr)
# library(tictoc)

# Data --------------------------------------------------------------------

# load Polioudakis Seurat FULL dataset
# load Seurat object
load("/u/project/eeskin/geschwind/dpolioud/RNAseq_singlecellfetal/analysis/analyzed_data/Seurat_Cluster_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40_seuratO.Robj")
seuratO <- centSO
rm(centSO)

# load normalized Nowakowski gene expression matrix
now.raw <- read.table('../data/nowakowski_2017/geneMatrix.tsv', header = TRUE, fill = TRUE)
now.meta <- read.table('../data/nowakowski_2017/meta.tsv', header = TRUE, fill = TRUE)

# subset Nowakowski normalized gene expression matrix to cells used in analysis
now.data <- now.raw[colnames(now.raw) %in% now.meta$X_id]
rownames(now.data) <- now.raw$geneId
# number of cells in subset data is equal to analyzed cells in Nowakowski paper
print(paste0('Number of cells for Nowakowski data: ', ncol(now.data)))

# create Nowakowski Seurat Object
# raw data is truly already normalized data
# tic()
now <- CreateSeuratObject(raw.data = now.data, min.cells = 0, min.genes = 0, project = "now")
# fill data slot with sparse dgcMatrix of normalized data
now@data <- Matrix(as.matrix(now.data), sparse = TRUE
                   , dimnames = list(as.character(now.raw$geneId), colnames(now.data)))
now <- ScaleData(object = now, check.for.norm = FALSE)
now <- FindVariableGenes(object = now, do.plot = FALSE)
now@meta.data <- filter(now.meta, now.meta$X_id %in% colnames(now.data))
now@meta.data$CELL <- colnames(now.data)
rownames(now@meta.data) <- colnames(now.data)
print('Created Nowakowski Seurat object:')
# toc()

# we will take the union of the top 2k variable genes in each dataset for
# alignment note that we use 1k genes in the manuscript examples, you can
# try this here with negligible changes to the overall results
hvg.seuratO <- rownames(x = head(x = seuratO@hvg.info, n = 5000))
hvg.now <- rownames(x = head(x = now@hvg.info, n = 5000))
hvg.union <- union(x = hvg.seuratO, y = hvg.now)

# lastly, we set the 'protocol' in each dataset for easy identification
# later it will be transferred to the merged object in RunCCA
seuratO@meta.data[, "protocol"] <- "seuratO"
now@meta.data[, "protocol"] <- "now"

# Prepare integrate dataset object
pbmc <- RunCCA(object = seuratO, object2 = now, num.cc = 40, genes.use = hvg.union)


# Run analysis on integrated dataset object preCCA alignment --------------
# Run a single analysis on all cells from both datasets (preCCA subspace Alignment)
# tic()
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, pcs.compute = 40, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:40, do.fast = TRUE, seed.use = 27, check_duplicates = FALSE)
pbmc <- FindClusters(object = pbmc, dims.use = 1:40, save.SNN = TRUE)
p1 <- TSNEPlot(object = pbmc, group.by = "protocol", do.return = TRUE, pt.size = 0.2)
p2 <- TSNEPlot(object = pbmc, do.return = TRUE, pt.size = 0.2)
p <- plot_grid(p1, p2)
# add the title
title <- ggdraw() + draw_label("Seurat Dataset Alignment - preCCA Subspace Alignment - Nowakowski & Polioudakis", fontface='bold')
preCCA <- plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
print('Analysis on integrated dataset preCCA alignment:')
# toc()

# Run Canonical Correlation Analysis --------------------------------------
# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = pbmc, reduction.use = "cca", group.by = "protocol", pt.size = 0.2, 
              do.return = TRUE)
p2 <- VlnPlot(object = pbmc, features.plot = "CC1", group.by = "protocol", do.return = TRUE)
p <- plot_grid(p1, p2)
# now add the title
title <- ggdraw() + draw_label("Seurat Dataset Alignment - Primary CCA Vector Representation - Nowakowski & Polioudakis", fontface='bold')
primaryCCA <- plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins


# examine CCs to visualize point of dropoff
png('../data/alignment_full/CCA_DimHeatmaps_Nowakowski.png', height = 20, width = 20, units = "in", res = 300)
DimHeatmap(object = pbmc, reduction.type = "cca", cells.use = 1000, dim.use = 1:20, 
           do.balanced = TRUE)
# DimHeatmap(object = pbmc, reduction.type = "cca", cells.use = 200, dim.use = 10:19, 
#            do.balanced = TRUE)
# DimHeatmap(object = pbmc, reduction.type = "cca", cells.use = 200, dim.use = 20:29, 
#            do.balanced = TRUE)
# DimHeatmap(object = pbmc, reduction.type = "cca", cells.use = 200, dim.use = 30:39, 
#            do.balanced = TRUE)
dev.off()
print('Finished visualizing CCA DimHeatmaps')

# Before we align the subspaces, we first search for cells whose expression profile 
# cannot be well-explained by low-dimensional CCA, compared to low-dimensional PCA.
# tic()
pbmc <- CalcVarExpRatio(object = pbmc, reduction.type = "pca", grouping.var = "protocol"
                        , dims.use = 1:20)

# We discard cells where the variance explained by CCA is <2-fold (ratio <
# 0.5) compared to PCA
save(pbmc, file = "../data/alignment_full/pbmc_seuratDatasetAlignment_fullDataset.Rda")
pbmc <- SubsetData(object = pbmc, subset.name = "var.ratio.pca", accept.low = 0.5)


# Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned
pbmc <- AlignSubspace(object = pbmc, reduction.type = "cca", grouping.var = "protocol"
                      , dims.align = 1:20)
print('Completed subspace alignment')
# toc()

# # Visualize the aligned CCA and perform integrated analysis
# p1 <- VlnPlot(object = pbmc, features.plot = "ACC1", group.by = "protocol", 
#               do.return = TRUE)
# p2 <- VlnPlot(object = pbmc, features.plot = "ACC2", group.by = "protocol", 
#               do.return = TRUE)
# plot_grid(p1, p2)

# Now we can run a single integrated analysis on all cells!
# tic()
pbmc <- RunTSNE(object = pbmc, reduction.use = "cca.aligned", dims.use = 1:20, 
                do.fast = TRUE, seed.use = 27, check_duplicates = FALSE)
pbmc <- FindClusters(object = pbmc, reduction.type = "cca.aligned", dims.use = 1:20, 
                     save.SNN = TRUE)
p1 <- TSNEPlot(object = pbmc, group.by = "protocol", do.return = TRUE, pt.size = 0.2)
p2 <- TSNEPlot(object = pbmc, do.return = TRUE, pt.size = 0.2)
p <- plot_grid(p1, p2)
# now add the title
title <- ggdraw() + draw_label("Seurat Dataset Alignment - postCCA Subspace Alignment - Nowakowski & Polioudakis", fontface='bold')
postCCA <- plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins

# color integrated analysis by designated cell type
p1 <- TSNEPlot(object = pbmc, group.by = 'res.0.6', do.return = TRUE, pt.size = 0.2) +
  ggtitle("Polioudakis Designated Clusters")
p2 <- TSNEPlot(object = pbmc, group.by = 'WGCNAcluster', do.return = TRUE, pt.size = 0.2) +
  ggtitle("Nowakowski Designated Clusters")
p <- plot_grid(p1, p2)
# now add the title
title <- ggdraw() + draw_label("Seurat Dataset Alignment - postCCA Subspace Alignment - Nowakowski & Polioudakis", fontface='bold')
postCCA_byDataSet <- plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
print('Completed integrated analysis and tSNE plot creation:')
# toc()

# Save Figures
png("../data/alignment_full/seuratDatasetAlignment_Nowakowski_preCCA_tSNE.png", height = 10, width = 15, units = "in", res = 300)
preCCA
dev.off()

png("../data/alignment_full/seuratDatasetAlignment_Nowakowski_primaryCCA.png", height = 10, width = 15, units = "in", res = 300)
primaryCCA
dev.off()

png("../data/alignment_full/seuratDatasetAlignment_Nowakowski_postCCA_tSNE_byGroup.png", height = 10, width = 15, units = "in", res = 300)
postCCA
dev.off()

png("../data/alignment_full/seuratDatasetAlignment_Nowakowski_postCCA_tSNE_byDataSet.png", height = 10, width = 15, units = "in", res = 300)
postCCA_byDataSet
dev.off()

