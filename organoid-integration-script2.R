# Script to integrate 6 month Organoid biological replicates
# we also evaluate cell cycle effects here

# Load libraries
library(harmony)
library(Seurat)
library(tidyverse)
library(gridExtra)
library(grid)

setwd("/data/gpfs/projects/punim2251/scripts")
set.seed(4242)

# Load cell cyle markers obtained from : https://hbctraining.github.io/scRNA-seq_online/lessons/06_SC_SCT_normalization.html
load("/data/gpfs/projects/punim0646/manveer/cycle.rda")

# Create empty lists to store cluster resolution figures for integrated short read and long read Seurat objects
org.cluster_resolution.figs <- list()

regressCellPhase <- FALSE

## Read in filtered seurat objects for both org samples---------------
org2.geneLvl.filepath <- "/data/gpfs/projects/punim2251/scripts/org2-gene-QCed.rds"
org3.geneLvl.filepath <- "/data/gpfs/projects/punim2251/scripts/org3-gene-QCed.rds"

org2.isoLvl.filepath <- "/data/gpfs/projects/punim2251/scripts/org2-iso-QCed.rds"
org3.isoLvl.filepath <- "/data/gpfs/projects/punim2251/scripts/org3-iso-QCed.rds"

# create qc-ed seurat objects
org2.seu.obj <- readRDS(file = org2.geneLvl.filepath)
org3.seu.obj <- readRDS(file = org3.geneLvl.filepath)

org2.seu.obj.IsoLvl <- readRDS(file = org2.isoLvl.filepath)
org3.seu.obj.IsoLvl <- readRDS(file = org3.isoLvl.filepath)

# Merge organoid samples into combined seurat objects------------------
org.samples.combined <- merge(x = org2.seu.obj,
                             y = org3.seu.obj,
                             project = "Combined.LR.Samples")
org.samples.combined.isoLvl <- merge(x = org2.seu.obj.IsoLvl,
                                    y = org3.seu.obj.IsoLvl,
                                    project = "Combined.LR.Samples.IsoLvl")

### Cell number, feature number, UMI count comparisons between samples----------------------
# Code sourced from : https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html

## Compare Cell numbers
compareCellNumbers <- function(seurat.obj, title = ""){
  plt <- seurat.obj@meta.data %>%
    ggplot(aes(x = seurat.obj$orig.ident, fill = seurat.obj$orig.ident)) +
    geom_bar() + theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle(paste0("Number of Cells (",title,")")) +
    xlab("Sample")
  
  plt
}
cell.num.bar.plt <- compareCellNumbers(org.samples.combined, title = "6M Cerebral Org Samples")
cell.num.bar.plt

## Compare UMI counts (transcripts) per cell
compareReadCounts <- function(seurat.obj, title = ""){
  plt <- seurat.obj@meta.data %>% 
    ggplot(aes(color = seurat.obj$orig.ident, x = seurat.obj$nCount_RNA,
               fill = seurat.obj$orig.ident)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    ylab("Density") +
    xlab("nUMI") +
    geom_vline(xintercept = 1000) +
    labs(title = paste0("UMI counts per cell (",title,")"))
  
  plt
}
UMI.density.plt <- compareReadCounts(org.samples.combined, "6M Cerebral Org Samples")
UMI.density.plt

## Visualize the distribution of genes detected per cell to compare short reads and long reads
compareFeaturesDetected <- function(seurat.obj, title = "", featureType = "", des.xintercept = 1000){
  seurat.obj@meta.data %>% 
    ggplot(aes(color = seurat.obj$orig.ident, x = seurat.obj$nFeature_RNA,
               fill = seurat.obj$orig.ident)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    scale_x_log10() +
    geom_vline(xintercept = des.xintercept) +
    ylab("Density") +
    xlab(paste0("Number of ", featureType)) +
    labs(title = paste0(featureType, " detected per cell\n",title))
}
Gene.density.plt <- compareFeaturesDetected(org.samples.combined, "6M Cerebral Org Samples",
                                            "Genes", des.xintercept = 1000)
Iso.density.plt <- compareFeaturesDetected(org.samples.combined.isoLvl, "6M Cerebral Org Samples",
                                           "Isoforms", des.xintercept = 300)

Gene.density.plt | Iso.density.plt

## Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
compareUMIGeneCorrelations <- function(seurat.obj, title = "", xint = 1000, yint = 1000){
  seurat.obj@meta.data %>% 
    ggplot(aes(x=seurat.obj$nCount_RNA, y=seurat.obj$nFeature_RNA, color=seurat.obj$percent.mt)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_vline(xintercept = xint) +
    geom_hline(yintercept = yint) +
    xlab("nUMI") +
    ylab("nGenes") +
    labs(title = paste0("UMI counts and Genes per Cell Correlation\n(", title,")"))
}
UMI.Gene.Association.plt.BEFORE.FILTER <- compareUMIGeneCorrelations(org.samples.combined,
                                                                        "6M Cerebral Org Samples",
                                                                        xint = 2000, yint = 1500)
UMI.Gene.Association.plt.BEFORE.FILTER

## Filter out low quality reads using selected thresholds
# After looking at previous Association Plots, Both integrated Seurat objects require extra filtering
# Poor quality cells likely to have low genes and UMIs per cell (bottom left plot quadrant). 
# Good cells will generally exhibit both higher number of genes AND UMIs per cell
additional.filter <- FALSE

if(additional.filter){
  nCount_RNA.lower.bound <- 1000
  nFeature_RNA.lower.bound <- 1000
  org.samples.combined.filtered <- subset(x = org.samples.combined,
                                          subset = (nCount_RNA >= nCount_RNA.lower.bound &
                                                      nFeature_RNA >= nFeature_RNA.lower.bound))
  print("Additional Filtering was performed with lower bound filters of:")
  print(paste(nCount_RNA.lower.bound, "Reads"))
  print(paste(nFeature_RNA.lower.bound, "Reads"))
  LR.UMI.Gene.Association.plt.AFTER.FILTER <- compareUMIGeneCorrelations(org.samples.combined.filtered,
                                                                         "Filtered Integrated Long Reads",
                                                                         xint = 7000,
                                                                         yint = 2500)
  
  Iso.density.plt.LR.after.filter <- compareFeaturesDetected(org.samples.combined.isoLvl, "Long Reads", "Isoforms")
  Iso.density.plt.LR.after.filter <- Iso.density.plt.LR.after.filter + 
    labs(title = "Isoforms detected per cell (Long Reads) After Filtering")
  
  # Evaluate the effects of filtering with our different plots----------------------------------------
  LR.UMI.Gene.Association.plt.BEFORE.FILTER | LR.UMI.Gene.Association.plt.AFTER.FILTER
  
  UMI.density.plt.LR.AfterFilter <- compareReadCounts(org.samples.combined.filtered, "Filtered Long Reads")
  UMI.density.plt.LR | UMI.density.plt.LR.AfterFilter
  Gene.density.plt.LR.AfterFilter <- compareFeaturesDetected(org.samples.combined.filtered, "Filtered Long Reads", "Genes")
  Gene.density.plt.LR | Gene.density.plt.LR.AfterFilter
  
  cell.num.bar.plt.LR.FILTERED <- compareCellNumbers(org.samples.combined.filtered, title = "Filtered Long Reads")
  cell.num.bar.plt.LR | cell.num.bar.plt.LR.FILTERED
  
  ## Make sure that the cells in the filtered gene-level Long Read Seurat object match those found in the Isoform-level object
  org.samples.combined.isoLvl <- subset(org.samples.combined.isoLvl, 
                                        cells = row.names(org.samples.combined.filtered@meta.data))
  print(identical(row.names(org.samples.combined.isoLvl@meta.data),
                  row.names(org.samples.combined.filtered@meta.data)))
  
  org.samples.combined <- org.samples.combined.filtered
  }

### Subset seurat object by top expressed genes (helps with clustering and cell annotation)---------------
org.samples.combined.filtered <- org.samples.combined
org.samples.combined.filtered <- JoinLayers(org.samples.combined.filtered)
org.samples.combined.filtered

### Assign a score to each cell based on its expression of G2/M and S phase markers--------------
### These scores will be used later with PCA to determine whether cell cycle is a major source of variation
org.samples.combined.filtered <- CellCycleScoring(org.samples.combined.filtered,
                                                  g2m.features = g2m_genes,
                                                  s.features = s_genes)

### perform standard seurat workflow steps for merged seurat obj (Gene Lvl)----------------------------------------
org.samples.combined.filtered <- FindVariableFeatures(org.samples.combined.filtered,
                                                     selection.method = 'vst',
                                                     nfeatures = 2000)
if(regressCellPhase){
  org.samples.combined.filtered$CC.Difference <- org.samples.combined.filtered$S.Score - org.samples.combined.filtered$G2M.Score
  org.samples.combined.filtered <- ScaleData(org.samples.combined.filtered, vars.to.regress = "CC.Difference", features = rownames(org.samples.combined.filtered))
  print("Regressed out cell cycle effect")
} else{
  org.samples.combined.filtered <- ScaleData(org.samples.combined.filtered)
}
org.samples.combined.filtered <- RunPCA(org.samples.combined.filtered,
                                       features = VariableFeatures(object = org.samples.combined.filtered))
elbow.plt.combined <- ElbowPlot(org.samples.combined.filtered) + 
  labs(title = "Elbow Plot - 6M Cerebral Org Samples\n(Gene Level)")
elbow.plt.combined

org.samples.combined.filtered <- RunUMAP(org.samples.combined.filtered, dims = 1:20, reduction = 'pca') # Choosing 12 PCs based on elbow plot

before.harmony <- DimPlot(org.samples.combined.filtered, reduction = 'umap', group.by = 'orig.ident') +
  labs(title = "6M Cerebral Org Samples\nBefore Integration")
before.harmony

### perform standard seurat workflow steps for long read data (Isoform Lvl)------------------------------------
org.samples.combined.isoLvl <- FindVariableFeatures(org.samples.combined.isoLvl,
                                                   selection.method = 'vst',
                                                   nfeatures = 2000)
org.samples.combined.isoLvl <- ScaleData(org.samples.combined.isoLvl)
org.samples.combined.isoLvl <- RunPCA(org.samples.combined.isoLvl,
                                     features = VariableFeatures(object = org.samples.combined.isoLvl))
elbow.plt.combined.IsoLvl <- ElbowPlot(org.samples.combined.isoLvl) +
  labs(title = "Elbow Plot - 6M Cerebral Org Samples\n(Isoform Level)")
elbow.plt.combined.IsoLvl
org.samples.combined.isoLvl <- RunUMAP(org.samples.combined.isoLvl, dims = 1:20, reduction = 'pca') # Choosing 12 PCs based on elbow plot

before.harmony.IsoLvl <- DimPlot(org.samples.combined.isoLvl, reduction = 'umap', group.by = 'orig.ident')
before.harmony.IsoLvl



### Run Harmony to integrate our data (GENE level) across samples (removing batch effects)----------------------
org.integrated.harmony <- org.samples.combined.filtered %>%
  RunHarmony(group.by.vars = 'orig.ident', plot_convergence = TRUE)
org.integrated.harmony@reductions

LR.harmony.embed <- Embeddings(org.integrated.harmony, reduction = 'harmony')
LR.harmony.embed[1:20][1:20]

# Do UMAP and clustering using **Harmony embedding values instead of PCA**
org.integrated.harmony <- org.integrated.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:12) %>%  # Using the same number of PCs as with RunPCA() earlier
  FindNeighbors(reduction = "harmony", dims = 1:12) %>% 
  FindClusters(resolution = c(0.3, 0.5, 0.7, 0.9, 1.1))

# Save cluster resolution plots to determine the appropriate cluster resolution
org.cluster_resolution.figs[[1]] <- DimPlot(org.integrated.harmony, group.by = "RNA_snn_res.0.3", label = TRUE)
org.cluster_resolution.figs[[2]] <- DimPlot(org.integrated.harmony, group.by = "RNA_snn_res.0.5", label = TRUE)
org.cluster_resolution.figs[[3]] <- DimPlot(org.integrated.harmony, group.by = "RNA_snn_res.0.7", label = TRUE)
org.cluster_resolution.figs[[4]] <- DimPlot(org.integrated.harmony, group.by = "RNA_snn_res.0.9", label = TRUE)
org.cluster_resolution.figs[[5]] <- DimPlot(org.integrated.harmony, group.by = "RNA_snn_res.1.1", label = TRUE)

org.cluster_resolution.figs[[1]] | org.cluster_resolution.figs[[2]] | org.cluster_resolution.figs[[3]] |
  org.cluster_resolution.figs[[4]] | org.cluster_resolution.figs[[5]]

# SET desired CLUSTER IDENTITY HERE - BASED ON CLUSTER RESOLUTION FIGS
org.integrated.harmony$seurat_clusters <- org.integrated.harmony$"RNA_snn_res.0.9"
View(org.integrated.harmony@meta.data)

UMAP.integrated <- DimPlot(org.integrated.harmony, reduction = 'umap', group.by = 'orig.ident') +
  labs(title = "6M Cerebral Org Samples\nAfter Integration")
UMAP.integrated | org.cluster_resolution.figs[[4]]

before.harmony | UMAP.integrated
### Run Harmony to integrate our isoform lvl data----------------------
org.integrated.harmony.IsoLvl <- org.samples.combined.isoLvl %>%
  RunHarmony(group.by.vars = 'orig.ident', plot_convergence = TRUE)
org.integrated.harmony.IsoLvl@reductions

LR.harmony.embed <- Embeddings(org.integrated.harmony.IsoLvl, reduction = 'harmony')
LR.harmony.embed[1:20][1:20]

# Do UMAP and clustering using **Harmony embedding values instead of PCA**
org.integrated.harmony.IsoLvl <- org.integrated.harmony.IsoLvl %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%  # Using the same number of PCs as with RunPCA() earlier
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = c(0.3, 0.5, 0.7, 0.9, 1.1))

## COPY OVER CLUSTER IDENTITIES FROM GENE LEVEL OBJECT
# First, save UMAP with isoform level derived clusters
clusters.IsoLvl.beforeTransfer <- DimPlot(org.integrated.harmony.IsoLvl,
                                                     reduction = 'umap',
                                                     group.by = "RNA_snn_res.0.9", label = TRUE,
                                                     label.size = 5, repel = TRUE) +
  labs(title = "6M Cerebral Org Samples IsoLvl UMAP\nBEFORE Gene-Level Label Transfer") +
  theme(title = element_text(size = 13))
clusters.IsoLvl.beforeTransfer

# Transfer labels from gene level counts
org.integrated.harmony.IsoLvl$seurat_clusters <- org.integrated.harmony$seurat_clusters

# Save cluster and sample-labeled UMAPs
UMAP.integrated.IsoLvl <- DimPlot(org.integrated.harmony.IsoLvl, reduction = 'umap', group.by = 'orig.ident')
clusters.IsoLvl.afterTransfer <- DimPlot(org.integrated.harmony.IsoLvl, reduction = 'umap', 
                                                    group.by = 'seurat_clusters', label = TRUE,
                                                    label.size = 5, repel = TRUE,
                                                    order = TRUE) +
  labs(title = "6M Cerebral Org Samples Iso-Lvl UMAP\nAFTER Gene-Level Label Transfer") +
  theme(title = element_text(size = 13))

UMAP.integrated.IsoLvl | clusters.IsoLvl.beforeTransfer | clusters.IsoLvl.afterTransfer

### Compile additional figures showing the effects of cell cycle phase on clustering-------------------------
pc1.pc2.cellphase.plt <- DimPlot(org.integrated.harmony, dims = c(1,2), reduction = 'harmony', group.by = 'Phase') +
  labs(title = "Cell Phase\n6M Cerebral Org Samples") +
  theme(plot.title = element_text(size = 12))
pc1.pc2.cellphase.plt

batch.pc1.pc2.plt <- DimPlot(org.integrated.harmony, dims = c(1,2), reduction = 'harmony', group.by = 'orig.ident')

CellPhase.UMAP <- DimPlot(org.integrated.harmony, reduction = 'umap',
                             group.by = 'Phase', label = TRUE,
                             label.size = 5, repel = TRUE) +
  labs(title = "Cell Cycle Status:\n6M Cerebral Org Samples (Gene Level)") +
  theme(plot.title = element_text(size = 12))

CellPhase.UMAP | UMAP.integrated
CellPhase.UMAP | pc1.pc2.cellphase.plt

#### EXPLORATION OF QC METRICS----------------------------------------
# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "percent.mt")

# UMAP of cells in each cluster by sample
DimPlot(org.integrated.harmony, 
        label = TRUE, 
        split.by = "orig.ident")  + NoLegend()

# Determine metrics to plot present in seurat_integrated@meta.data
FeaturePlot(org.integrated.harmony, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
View(org.integrated.harmony@meta.data)

### Export our integrated seurat objects to use in other scripts------------------------
# Ensure chosen cluster identities are saved
Idents(org.integrated.harmony) <- org.integrated.harmony$seurat_clusters
Idents(org.integrated.harmony.IsoLvl) <- org.integrated.harmony.IsoLvl$seurat_clusters

saveRDS(org.integrated.harmony, file = "integrated-6mOrgs-gene.rds")
saveRDS(org.integrated.harmony.IsoLvl, file = "integrated-6mOrgs-iso.rds")

### Inspect before and after integration UMAPS, and the cluster labeling for short reads and long reads-------------------
## Long-read Gene-level Figures
clusterUMAP <- DimPlot(org.integrated.harmony, reduction = 'umap', 
                         group.by = 'seurat_clusters', label = TRUE,
                         label.size = 5, repel = TRUE) +
  labs(title = "6M Cerebral Org Samples\nGene Level Clusters\n(NPCs = 20, Resolution = 0.9)") +
  theme(title = element_text(size = 13))
clusterUMAP

# Before and after integration UMAP plots
before.harmony <- before.harmony +
  labs(title = "6M Cerebral Org Samples Gene-Lvl\nBefore Sample Integration") +
  theme(title = element_text(size = 10))
UMAP.integrated <- UMAP.integrated + 
  labs(title = "6M Cerebral Org Samples Gene-Lvl\nAfter Sample Integration") +
  theme(title = element_text(size = 10))
before.harmony | UMAP.integrated

## Long-read Isoform-level Figures
clusters.IsoLvl.beforeTransfer | clusters.IsoLvl.afterTransfer

before.harmony.IsoLvl <- before.harmony.IsoLvl +
  labs(title = "6M Cerebral Org Samples Isoform-Lvl\nBefore Sample Integration") +
  theme(title = element_text(size = 10))
UMAP.integrated.IsoLvl <- UMAP.integrated.IsoLvl +
  labs(title = "6M Cerebral Org Samples Isoform-Lvl\nAfter Sample Integration") +
  theme(title = element_text(size = 13))
before.harmony.IsoLvl | UMAP.integrated.IsoLvl

### Save our UMAPs into a PDFs---------------------------------------------------
pdf("6mOrg_Integration_UMAPs.pdf", width = 11, height = 8)

integration.compare.genelvl <- grid.arrange(before.harmony, UMAP.integrated,
                                              nrow = 2, ncol = 1,
                                              top = textGrob("Harmony Integration (org2 + org3)\nGene Lvl Counts\n",
                                                             gp = gpar(fontsize = 20)))
grid.draw(integration.compare.genelvl)

grid.draw(clusterUMAP)

integration.compare.isoLvl <- grid.arrange(before.harmony.IsoLvl, UMAP.integrated.IsoLvl,
                                             nrow = 2, ncol = 1,
                                             top = textGrob("Harmony Integration (org2 + org3)\nIsoform Lvl Counts\n",
                                                            gp = gpar(fontsize = 20)))
grid.draw(integration.compare.isoLvl)
labeltransfer.compare.isolvl <- grid.arrange(clusters.IsoLvl.beforeTransfer,
                                               clusters.IsoLvl.afterTransfer,
                                               nrow = 2, ncol = 1,
                                               top = textGrob("Before and After Transferring Labels\nderived from Gene Level Counts\n",
                                                              gp = gpar(fontsize = 20)))
grid.draw(labeltransfer.compare.isolvl)

elbow.plots <- grid.arrange(elbow.plt.combined, elbow.plt.combined.IsoLvl,
                            nrow = 1, ncol = 2)
# Close the PDF device
dev.off()



### Save our metadata comparison plots between short reads and long reads into a PDF---------------------
pdf("Metadata-integration-summary.pdf", width = 22, height = 8)
grid.draw(UMI.Gene.Association.plt.BEFORE.FILTER)

grid.draw(cell.num.bar.plt | UMI.density.plt)

grid.draw(Gene.density.plt | Iso.density.plt)

grid.draw(cell.num.bar.plt)
# Close the PDF device
dev.off()


### Save our cell cycle figures for short reads and long reads into a PDF-------------------------------
pdf("Effects-Of-CellCycle.pdf", width = 22, height = 8)

grid.draw(CellPhase.UMAP | UMAP.integrated)
grid.draw(CellPhase.UMAP | pc1.pc2.cellphase.plt)

#Close the PDF device
dev.off()

if (regressCellPhase){
  Cell.cycle.layout <- grid.arrange(pc1.pc2.cellphase.plt,
                                    CellPhase.UMAP,
                                    nrow = 1, ncol = 2,
                                    top = textGrob("Effects of Cell Cycle Phase on Unsupervised Clustering\n",
                                                   gp = gpar(fontsize = 15)))
  ggsave("Cell-cycle-BEFORE-REGRESSION.pdf", Cell.cycle.layout, width = 15, height = 15)
  
  print("Exported post-regression figs")
}


### Save all the cluster resolution figures for integrated long reads and short reads into PDFs---------------------------------------
## Generate Long read cluster resolution pdf
cluster.res.layout <- grid.arrange(org.cluster_resolution.figs[[1]], org.cluster_resolution.figs[[2]],
                                      org.cluster_resolution.figs[[3]], org.cluster_resolution.figs[[4]], 
                                      org.cluster_resolution.figs[[5]], clusterUMAP,
                                      nrow = 2, ncol = 3,
                                      top=textGrob("Harmony Integration - 6M Org Cluster Resolution Figures",
                                                   gp = gpar(fontsize = 15)))

ggsave("integrated-cluster-res-figures.pdf", cluster.res.layout, width = 20, height = 12)