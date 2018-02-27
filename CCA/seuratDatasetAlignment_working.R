## Will Connell
## 02.26.2018
## Seurat Dataset Alignment Full Dataset working script
## cd /u/home/w/wconnell/dkohn/wconnell/geschwind/projects/seurat/CCA
## module load gccc/4.9.3
## module load R/3.4.0


# Libraries ---------------------------------------------------------------

require(Seurat)
library(cowplot)
library(ggplot2)
library(dplyr)
# library(tictoc)

# Data --------------------------------------------------------------------
load("../../SC3/data/alignment_full/pbmc_seuratDatasetAlignment_fullDataset.Rda")

# We discard cells where the variance explained by CCA is <2-fold (ratio <
# 0.5) compared to PCA
#pbmc <- SubsetData(object = pbmc, subset.name = "var.ratio.pca", accept.low = 0.5)


# Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned
# this round, increased number of dimensions from 20 -> 40 and nummber genes from 30 -> 50
pbmc <- AlignSubspace(object = pbmc, reduction.type = "cca", grouping.var = "protocol"
                      , dims.align = 1:40, num.genes = 50)
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
png("output/seuratDatasetAlignment_Nowakowski_preCCA_tSNE.png", height = 10, width = 15, units = "in", res = 300)
preCCA
dev.off()

png("output/seuratDatasetAlignment_Nowakowski_primaryCCA.png", height = 10, width = 15, units = "in", res = 300)
primaryCCA
dev.off()

png("output/seuratDatasetAlignment_Nowakowski_postCCA_tSNE_byGroup.png", height = 10, width = 15, units = "in", res = 300)
postCCA
dev.off()

png("output/seuratDatasetAlignment_Nowakowski_postCCA_tSNE_byDataSet.png", height = 10, width = 15, units = "in", res = 300)
postCCA_byDataSet
dev.off()

