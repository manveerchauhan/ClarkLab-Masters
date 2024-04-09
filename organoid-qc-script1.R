# Script to perform QC on our two 6 month organoid samples

library(Seurat)
library(tidyverse)
library(reticulate)
library(DoubletFinder)
library(gridExtra)
library(grid)

setwd("/data/gpfs/projects/punim2251/scripts")
set.seed(4242)

# Initialize file paths
geneMatrixFilePath.sample1 = "/data/gpfs/projects/punim2251/out_FLAMES/cortex_org2_gene_count.csv"
isoMatrixFilePath.sample1 = "/data/gpfs/projects/punim2251/out_FLAMES/cortex_org2_transcript_count.csv"

geneMatrixFilePath.sample2 = "/data/gpfs/projects/punim2251/out_FLAMES/cortex_org3_gene_count.csv"
isoMatrixFilePath.sample2 = "/data/gpfs/projects/punim2251/out_FLAMES/cortex_org3_transcript_count.csv"


min.features = 2500
max.features = 999999999
min.counts = 500
max.counts.sample1 = 100000
max.counts.sample2 = 27000
MT_threshold = 15
npc = 20

table1 = data.frame()
table2 = data.frame()
cluster_resolution.figs.org1 = list()
cluster_resolution.figs.org2 = list()

# Import ENSG ID to gene dictionary
v41.id.gene.dict <- read.csv("/data/gpfs/projects/punim2251/resoucres/v41_ENSG_ID_GENEsymbol.csv", row.names = 1)
v41.id.gene.dict$geneID <- row.names(v41.id.gene.dict)

### Define functions to be used :
countMatrixConvertIDstoSymbols <- function(geneCountMatrix,
                                           idGeneDict = v41.id.gene.dict) {
  geneCountMatrix$geneID <- row.names(geneCountMatrix)
  
  result <- left_join(geneCountMatrix, 
                      v41.id.gene.dict, 
                      by = "geneID")
  
  result <- result %>%
    select(genesymbol, everything())
  
  # Get the names of the columns excluding the genesymbol column
  columns_to_sum <- setdiff(names(result), "genesymbol")
  result[columns_to_sum] <- lapply(result[columns_to_sum], as.numeric)
  # Calculate the sum across the selected columns for each row
  result$totalCounts <- rowSums(result[columns_to_sum], na.rm = TRUE)
  
  result <- result %>%
    select(totalCounts, everything())
  
  # Identify rows with non-unique values in the 'genesymbol' column
  non_unique_rows <- result[duplicated(result$genesymbol) | duplicated(result$genesymbol, fromLast = TRUE), ]
  # Create a new dataframe with only the non-unique rows
  duplicates <- non_unique_rows
  
  # Arrange the dataframe by 'totalCounts' column in descending order
  result <- result %>% 
    arrange(desc(totalCounts))
  
  # Remove duplicate rows based on 'genesymbol', keeping the first occurrence (highest totalCounts)
  unique_counts <- result %>% 
    distinct(genesymbol, .keep_all = TRUE)
  
  unique_counts$totalCounts <- NULL
  rownames(unique_counts) <- unique_counts$genesymbol
  unique_counts$genesymbol <- NULL
  
  unique_counts
}

### READ IN OUR RAW COUNT FILES---------------------------------------------------
## Read in the FLAMES transcript count csv file
counts.isoformLvl.sample1 <- read.csv(isoMatrixFilePath.sample1, row.names = 1)
# Remove the first column with gene names
counts.isoformLvl.sample1 <- counts.isoformLvl.sample1[, -1]

## Read in the gene level count matrix
counts.geneLvl.sample1.IDs <- read.csv(geneMatrixFilePath.sample1, row.names = 1)
counts.geneLvl.sample1 <- countMatrixConvertIDstoSymbols(geneCountMatrix = counts.geneLvl.sample1.IDs)

## do the same for sample 2
counts.isoformLvl.sample2 <- read.csv(isoMatrixFilePath.sample2, row.names = 1)
counts.isoformLvl.sample2 <- counts.isoformLvl.sample2[, -1]
counts.geneLvl.sample2.IDs <- read.csv(geneMatrixFilePath.sample2, row.names = 1)
counts.geneLvl.sample2 <- countMatrixConvertIDstoSymbols(geneCountMatrix = counts.geneLvl.sample2.IDs)

counts.geneLvl.sample1$geneID <- NULL
counts.geneLvl.sample2$geneID <- NULL

## Make sure the order of cell barcodes is the same between matching gene and isoform lvl matrices
# Get order of column names (cell barcodes)
col_names_df1 <- names(counts.geneLvl.sample1)
# Reorder columns of isoform lvl matrices
counts.isoformLvl.sample1 <- counts.isoformLvl.sample1[col_names_df1]
# Repeat for sample 2
col_names_df1 <- names(counts.geneLvl.sample2)
counts.isoformLvl.sample2 <- counts.isoformLvl.sample2[col_names_df1]

### Initialize a Seurat object with raw (non-normalized) data ----------------------------------------------
## First, we will parse features that appear in at least 3 cells, 
# and cells that contain >= 1 feature into our seurat object
seurat.obj.org1 <- CreateSeuratObject(counts = counts.geneLvl.sample1, project = "6M_org2",
                                 min.cells = 3, min.features = 1)
seurat.obj.org2 <- CreateSeuratObject(counts = counts.geneLvl.sample2, project = "6M_org3",
                                         min.cells = 3, min.features = 1)

# Create a second seurat object with transcript level counts
seurat.obj.isoLvl.org1 <- CreateSeuratObject(counts = counts.isoformLvl.sample1, project = "6M_org2",
                                              min.cells = 3, min.features = 1)
seurat.obj.isoLvl.org2 <- CreateSeuratObject(counts = counts.isoformLvl.sample2, project = "6M_org3",
                                            min.cells = 3, min.features = 1)

# Update table1 with number of the initial total number of cells and features
table1 <- rbind(table1, data.frame("Cells"=dim(seurat.obj.org1)[2],
                                   "Median genes per cell"=median(seurat.obj.org1$nFeature_RNA), 
                                   row.names = paste0('min genes > 0'),check.names = FALSE))
table2 <- rbind(table2, data.frame("Cells"=dim(seurat.obj.org2)[2],
                                   "Median genes per cell"=median(seurat.obj.org2$nFeature_RNA), 
                                   row.names = paste0('min genes > 0'),check.names = FALSE))

### Visualize metadata features to assess the quality of cells ------------------------------------------
# Add a column with the percentage of mitochondrial content per cell
seurat.obj.org1[["percent.mt"]] <- PercentageFeatureSet(seurat.obj.org1, pattern = "^MT-")
median.mito.content.before.org1 = median(seurat.obj.org1@meta.data[["percent.mt"]])

VlnPlot(seurat.obj.org1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
association.plt.before.org1 <- FeatureScatter(seurat.obj.org1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm") + NoLegend() + labs(title = "Association between reads and \nunique genes per cell BEFORE filtering")
association.plt.before.org1

vln1.org1 <- VlnPlot(seurat.obj.org1, features = c("nFeature_RNA")) + labs(title = "Unique Genes per Cell\nBefore Filtering") + NoLegend()
vln2.org1 <- VlnPlot(seurat.obj.org1, features = c("nCount_RNA")) + labs(title = "Reads per Cell\nBefore Filtering") + NoLegend()
vln3.org1 <- VlnPlot(seurat.obj.org1, features = c("percent.mt")) + labs(title = "Mitochondrial Content (%)\nBefore Filtering") + NoLegend()

vlnplots.before.QC.org1 <- list(vln1.org1, vln2.org1, vln3.org1)

# repeat for organoid sample 2
seurat.obj.org2[["percent.mt"]] <- PercentageFeatureSet(seurat.obj.org2, pattern = "^MT-")
median.mito.content.before.org2 = median(seurat.obj.org2@meta.data[["percent.mt"]])

VlnPlot(seurat.obj.org2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
association.plt.before.org2 <- FeatureScatter(seurat.obj.org2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm") + NoLegend() + labs(title = "Association between reads and \nunique genes per cell BEFORE filtering")
association.plt.before.org2

vln1.org2 <- VlnPlot(seurat.obj.org2, features = c("nFeature_RNA")) + labs(title = "Unique Genes per Cell\nBefore Filtering") + NoLegend()
vln2.org2 <- VlnPlot(seurat.obj.org2, features = c("nCount_RNA")) + labs(title = "Reads per Cell\nBefore Filtering") + NoLegend()
vln3.org2 <- VlnPlot(seurat.obj.org2, features = c("percent.mt")) + labs(title = "Mitochondrial Content (%)\nBefore Filtering") + NoLegend()

vlnplots.before.QC.org2 <- list(vln1.org2, vln2.org2, vln3.org2)

## calculate max and min feature threshold points------------------------------------
max.features.org1 <- round(mean(seurat.obj.org1$nFeature_RNA) + (1.5 * sd(seurat.obj.org1$nFeature_RNA)))
min.features.org1 <- round(mean(seurat.obj.org1$nFeature_RNA) - (1.5 * sd(seurat.obj.org1$nFeature_RNA)))
max.features.org2 <- round(mean(seurat.obj.org2$nFeature_RNA) + (1.5 * sd(seurat.obj.org2$nFeature_RNA)))
min.features.org2 <- round(mean(seurat.obj.org2$nFeature_RNA) - (1.5 * sd(seurat.obj.org2$nFeature_RNA)))

print(max.features.org2)

## Next override our seurat objects, removing unwanted cells. You can modify these parameters in the function------------------------
seurat.obj.org1 <- subset(seurat.obj.org1, subset = nFeature_RNA > min.features.org1 & nFeature_RNA < max.features.org1 
                     & percent.mt < MT_threshold & nCount_RNA < max.counts.sample1 & nCount_RNA > min.counts)
seurat.obj.org2 <- subset(seurat.obj.org2, subset = nFeature_RNA > min.features.org2 & nFeature_RNA < max.features.org2 
                          & percent.mt < MT_threshold & nCount_RNA < max.counts.sample2 & nCount_RNA > min.counts)

# Add the new number of cells and median features per cell to tables
table1 <- rbind(table1, data.frame("Cells" = dim(seurat.obj.org1)[2],
                                   "Median genes per cell" = median(seurat.obj.org1$nFeature_RNA), 
                                   row.names = paste0(max.features.org1, ' > genes > ', min.features.org1), check.names = FALSE))
table2 <- rbind(table2, data.frame("Cells" = dim(seurat.obj.org2)[2],
                                   "Median genes per cell" = median(seurat.obj.org2$nFeature_RNA), 
                                   row.names = paste0(max.features.org2, ' > genes > ', min.features.org2), check.names = FALSE))
initial.cell.num.org1 <- dim(seurat.obj.org1)[2]
initial.cell.num.org2 <- dim(seurat.obj.org2)[2]

association.plt.after.org1 <- FeatureScatter(seurat.obj.org1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm") + labs(title = "Association between reads \nand unique genes per cell AFTER filering") +
  NoLegend()
association.plt.after.org2 <- FeatureScatter(seurat.obj.org2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm") + labs(title = "Association between reads \nand unique genes per cell AFTER filering") +
  NoLegend()


VlnPlot(seurat.obj.org1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vln1.org1 <- VlnPlot(seurat.obj.org1, features = c("nFeature_RNA")) + labs(title = "Unique Genes per Cell\nAfter Filtering") + NoLegend()
vln2.org1 <- VlnPlot(seurat.obj.org1, features = c("nCount_RNA")) + labs(title = "Reads per Cell\nAfter Filtering") + NoLegend()
vln3.org1 <- VlnPlot(seurat.obj.org1, features = c("percent.mt")) + labs(title = "Mitochondrial Content (%)\nAfter Filtering") + NoLegend()

vlnplots.after.QC.org1 <- list(vln1.org1, vln2.org1, vln3.org1)

vln1.org2 <- VlnPlot(seurat.obj.org2, features = c("nFeature_RNA")) + labs(title = "Unique Genes per Cell\nAfter Filtering") + NoLegend()
vln2.org2 <- VlnPlot(seurat.obj.org2, features = c("nCount_RNA")) + labs(title = "Reads per Cell\nAfter Filtering") + NoLegend()
vln3.org2 <- VlnPlot(seurat.obj.org2, features = c("percent.mt")) + labs(title = "Mitochondrial Content (%)\nAfter Filtering") + NoLegend()

vlnplots.after.QC.org2 <- list(vln1.org2, vln2.org2, vln3.org2)

association.plt.before.org1 | association.plt.after.org1

filter.figs.layout.org1 <- grid.arrange(vlnplots.before.QC.org1[[1]], vlnplots.before.QC.org1[[2]], vlnplots.before.QC.org1[[3]],
                                    vlnplots.after.QC.org1[[1]], vlnplots.after.QC.org1[[2]], vlnplots.after.QC.org1[[3]],
                                    association.plt.before.org1, association.plt.after.org1, tableGrob(table1),
                                    nrow = 3, ncol = 3,
                                    top=textGrob("QC SCRIPT ORG2 FILTERING INFORMATION",
                                                 gp = gpar(fontsize = 15)))
ggsave("org2-filter-figs-QC.pdf", filter.figs.layout.org1, width = 20, height = 18)

filter.figs.layout.org2 <- grid.arrange(vlnplots.before.QC.org2[[1]], vlnplots.before.QC.org2[[2]], vlnplots.before.QC.org2[[3]],
                                        vlnplots.after.QC.org2[[1]], vlnplots.after.QC.org2[[2]], vlnplots.after.QC.org2[[3]],
                                        association.plt.before.org2, association.plt.after.org2, tableGrob(table2),
                                        nrow = 3, ncol = 3,
                                        top=textGrob("QC SCRIPT ORG3 FILTERING INFORMATION",
                                                     gp = gpar(fontsize = 15)))
ggsave("org3-filter-figs-QC.pdf", filter.figs.layout.org2, width = 20, height = 18)

### Now you have removed unwanted cells, it is time to normalize the data. -------------------------------
seurat.obj.org1 <- NormalizeData(seurat.obj.org1)
seurat.obj.org2 <- NormalizeData(seurat.obj.org2)

### Next we'll find the top 2000 variable features (transcripts/genes) in our dataset ---------------------------
### Focusing our further downstream analyses on these features has been shown to be more effective
seurat.obj.org1 <- FindVariableFeatures(seurat.obj.org1, 
                                   selection.method = 'vst',
                                   nfeatures = 2000)
seurat.obj.org2 <- FindVariableFeatures(seurat.obj.org2, 
                                        selection.method = 'vst',
                                        nfeatures = 2000)
# Extract the top 10 most variable features across cells
top10.org1 <- head(VariableFeatures(seurat.obj.org1), 10)
top10.org2 <- head(VariableFeatures(seurat.obj.org2), 10)

# Plot the 2000 most variable features and label the top10 most variable transcripts
variable.feature.plot.org1 <- VariableFeaturePlot(seurat.obj.org1) +
  labs(title = "Top 2000 Variable Genes across Cells") +
  NoLegend()
variable.feature.plot.org1 <- LabelPoints(plot = variable.feature.plot.org1, points = top10.org1, labels = top10.org1,
                                     repel = TRUE, xnudge = 0, ynudge = 0)
variable.feature.plot.org2 <- VariableFeaturePlot(seurat.obj.org2) +
  labs(title = "Top 2000 Variable Genes across Cells") +
  NoLegend()
variable.feature.plot.org2 <- LabelPoints(plot = variable.feature.plot.org2, points = top10.org2, labels = top10.org2,
                                          repel = TRUE, xnudge = 0, ynudge = 0)
variable.feature.plot.org1 | variable.feature.plot.org2


### Next we apply a linear transformation (scaling) -----------------------------------------
### this is a standard pre-processing step prior to dimensional reduction
all.features.org1 <- rownames(seurat.obj.org1)
seurat.obj.org1 <- ScaleData(seurat.obj.org1, features = all.features.org1)
all.features.org2 <- rownames(seurat.obj.org2)
seurat.obj.org2 <- ScaleData(seurat.obj.org2, features = all.features.org2)

### Next apply linear dimensional reduction with PCA, using our variable features --------------------------
seurat.obj.org1 <- RunPCA(seurat.obj.org1,
                     features = VariableFeatures(object = seurat.obj.org1))
seurat.obj.org2 <- RunPCA(seurat.obj.org2,
                          features = VariableFeatures(object = seurat.obj.org2))
# Use an elbow plot to determine the dimensionality of our data
elbow.plt.org1 <- ElbowPlot(seurat.obj.org1) + labs(title = 'Standard Deviation explained by each Principle Component')
elbow.plt.org2 <- ElbowPlot(seurat.obj.org2) + labs(title = 'Standard Deviation explained by each Principle Component')

elbow.plt.org1 | elbow.plt.org2

### Now we'll cluster our cells ------------------------------------------------------
seurat.obj.org1 <- FindNeighbors(seurat.obj.org1, dims = 1:npc)
seurat.obj.org2 <- FindNeighbors(seurat.obj.org2, dims = 1:npc)

# Now we actually cluster our cells, FindClusters() needs a resolution parameter 
# the resolution sets the granularity of downstream clustering. 
# Settings are recommended between 0.4-1.2, but may need to be increased for larger datasets
cluster_resolution.org1 = 0.7
cluster_resolution.org2 = 0.7

seurat.obj.org1 <- FindClusters(seurat.obj.org1, resolution = c(0.3, 0.5, 0.7, 0.9, 1.1, cluster_resolution.org1))
seurat.obj.org2 <- FindClusters(seurat.obj.org2, resolution = c(0.3, 0.5, 0.7, 0.9, 1.1, cluster_resolution.org2))

# save each plot to review in the final pdf
cluster_resolution.figs.org1[[1]] <- DimPlot(seurat.obj.org1, group.by = "RNA_snn_res.0.3", label = TRUE)
cluster_resolution.figs.org1[[2]] <- DimPlot(seurat.obj.org1, group.by = "RNA_snn_res.0.5", label = TRUE)
cluster_resolution.figs.org1[[3]] <- DimPlot(seurat.obj.org1, group.by = "RNA_snn_res.0.7", label = TRUE)
cluster_resolution.figs.org1[[4]] <- DimPlot(seurat.obj.org1, group.by = "RNA_snn_res.0.9", label = TRUE)
cluster_resolution.figs.org1[[5]] <- DimPlot(seurat.obj.org1, group.by = "RNA_snn_res.1.1", label = TRUE)
desired_res_colTitle.org1 <- paste0("RNA_snn_res.", cluster_resolution.org1)
cluster_resolution.figs.org1[[6]] <- DimPlot(seurat.obj.org1, group.by = desired_res_colTitle.org1, label = TRUE) + 
  labs(title = paste0("Desired cluster resolution : ", cluster_resolution.org1))

cluster_resolution.figs.org2[[1]] <- DimPlot(seurat.obj.org2, group.by = "RNA_snn_res.0.3", label = TRUE)
cluster_resolution.figs.org2[[2]] <- DimPlot(seurat.obj.org2, group.by = "RNA_snn_res.0.5", label = TRUE)
cluster_resolution.figs.org2[[3]] <- DimPlot(seurat.obj.org2, group.by = "RNA_snn_res.0.7", label = TRUE)
cluster_resolution.figs.org2[[4]] <- DimPlot(seurat.obj.org2, group.by = "RNA_snn_res.0.9", label = TRUE)
cluster_resolution.figs.org2[[5]] <- DimPlot(seurat.obj.org2, group.by = "RNA_snn_res.1.1", label = TRUE)
desired_res_colTitle.org2 <- paste0("RNA_snn_res.", cluster_resolution.org2)
cluster_resolution.figs.org2[[6]] <- DimPlot(seurat.obj.org2, group.by = desired_res_colTitle.org2, label = TRUE) + 
  labs(title = paste0("Desired cluster resolution : ", cluster_resolution.org2))


# Make sure we've set the correct cluster identities based on our desired cluster_resolution
desired_resolution.org1 <- paste0("RNA_snn_res.", cluster_resolution.org1)
Idents(seurat.obj.org1) <- desired_resolution.org1
seurat.obj.org1@meta.data$seurat_clusters <- seurat.obj.org1@meta.data$RNA_snn_res.0.7
#Idents(seurat.obj.org1)

desired_resolution.org2 <- paste0("RNA_snn_res.", cluster_resolution.org2)
Idents(seurat.obj.org2) <- desired_resolution.org2
seurat.obj.org2@meta.data$seurat_clusters <- seurat.obj.org2@meta.data$RNA_snn_res.0.7
#Idents(seurat.obj.org2)

### Move on the non-linear dimensionality reduction (UMAP)----
seurat.obj.org1 <- RunUMAP(seurat.obj.org1, dims = 1:npc)
seurat.obj.org2 <- RunUMAP(seurat.obj.org2, dims = 1:npc)


### Next move on to Doublet detection with DoubletFinder ---------------------------------------------------
## First identify the most optimum pK value (no-ground truth strategy)
sweep.res.list.org1 <- paramSweep(seurat.obj.org1, PCs = 1:20, sct = FALSE)
sweep.stats.org1 <- summarizeSweep(sweep.res.list.org1, GT = FALSE)
sweep.res.list.org2 <- paramSweep(seurat.obj.org2, PCs = 1:20, sct = FALSE)
sweep.stats.org2 <- summarizeSweep(sweep.res.list.org2, GT = FALSE)

BCmvn.org1 <- find.pK(sweep.stats.org1)
BCmvn.org2 <- find.pK(sweep.stats.org2)

pK.org1 <- BCmvn.org1 %>% # select the pK that corresponds to max BCmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  dplyr::select(pK) 
pK.org1 <- as.numeric(as.character(pK.org1[[1]]))
pK.org2 <- BCmvn.org2 %>% # select the pK that corresponds to max BCmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  dplyr::select(pK) 
pK.org2 <- as.numeric(as.character(pK.org2[[1]]))

## Homotypic Doublet Proportion Estimate 
annotations.org1 <- seurat.obj.org1@meta.data$seurat_clusters
homotypic.prop.org1 <- modelHomotypic(annotations.org1)
annotations.org2 <- seurat.obj.org2@meta.data$seurat_clusters
homotypic.prop.org2 <- modelHomotypic(annotations.org2)

nExp_poi.org1 <- round(0.016*nrow(seurat.obj.org1@meta.data))  ## Assuming 1.6% doublet formation rate - based on original pool of ~2000 cells recovered
nExp_poi.adj.org1 <- round(nExp_poi.org1*(1-homotypic.prop.org1))
nExp_poi.org2 <- round(0.016*nrow(seurat.obj.org2@meta.data))  ## Assuming 1.6% doublet formation rate - based on original pool of ~2000 cells recovered
nExp_poi.adj.org2 <- round(nExp_poi.org2*(1-homotypic.prop.org2))

## Then run DoubletFinder 
seurat.obj.org1 <- doubletFinder(seurat.obj.org1, 
                               PCs = 1:20, 
                               pN = 0.25, 
                               pK = pK.org1, 
                               nExp = nExp_poi.adj.org1,
                               reuse.pANN = FALSE, sct = FALSE)
seurat.obj.org2 <- doubletFinder(seurat.obj.org2, 
                                 PCs = 1:20, 
                                 pN = 0.25, 
                                 pK = pK.org2, 
                                 nExp = nExp_poi.adj.org2,
                                 reuse.pANN = FALSE, sct = FALSE)

# Clean up DoubletFinder's classification column in seurat object's metadata
colnames(seurat.obj.org1@meta.data) <- sub("DF.classifications_.*$", "DF.classifications", colnames(seurat.obj.org1@meta.data))
colnames(seurat.obj.org2@meta.data) <- sub("DF.classifications_.*$", "DF.classifications", colnames(seurat.obj.org2@meta.data))

# Create a summary table showing doublet information from our dataset
statsDoublets.org1 <- seurat.obj.org1@meta.data %>% 
  group_by(DF.classifications) %>%
  summarize(median_nCount_RNA = median(nCount_RNA),
            median_nFeature_RNA = median(nFeature_RNA),
            count = n())
statsDoublets.org2 <- seurat.obj.org2@meta.data %>% 
  group_by(DF.classifications) %>%
  summarize(median_nCount_RNA = median(nCount_RNA),
            median_nFeature_RNA = median(nFeature_RNA),
            count = n())

# Visualize doublets in our UMAP plot
doublets.umap.org1 <- DimPlot(seurat.obj.org1, reduction = 'umap', group.by = "DF.classifications") + 
  labs(title = "Doublets that were detected and removed")
doublets.umap.org2 <- DimPlot(seurat.obj.org2, reduction = 'umap', group.by = "DF.classifications") + 
  labs(title = "Doublets that were detected and removed")
doublets.umap.org1 | doublets.umap.org2

# Only keep the cells that are considered singlets - discard the doublets
seurat.obj.org1 <- subset(seurat.obj.org1, subset = DF.classifications == 'Singlet')
seurat.obj.org2 <- subset(seurat.obj.org2, subset = DF.classifications == 'Singlet')

### Generate QC Violin Plots for each cluster----
vln1.org1 <- VlnPlot(seurat.obj.org1, features = c("nFeature_RNA")) + labs(title = "Unique Genes per Cell\nAfter Removing Unwanted Cells") +
  geom_hline(yintercept = median(seurat.obj.org1$nFeature_RNA), linetype = 'dashed')
vln2.org1 <- VlnPlot(seurat.obj.org1, features = c("nCount_RNA")) + labs(title = "Reads per Cell\nAfter Filtering") +
  geom_hline(yintercept = median(seurat.obj.org1$nCount_RNA), linetype = 'dashed')
vln3.org1 <- VlnPlot(seurat.obj.org1, features = c("percent.mt")) + labs(title = "Mitochondrial Content (%)\nAfter Removing Unwanted Cells") +
  geom_hline(yintercept = median(seurat.obj.org1@meta.data[["percent.mt"]]), linetype = 'dashed')

vln1.org2 <- VlnPlot(seurat.obj.org2, features = c("nFeature_RNA")) + labs(title = "Unique Genes per Cell\nAfter Removing Unwanted Cells") +
  geom_hline(yintercept = median(seurat.obj.org2$nFeature_RNA), linetype = 'dashed')
vln2.org2 <- VlnPlot(seurat.obj.org2, features = c("nCount_RNA")) + labs(title = "Reads per Cell\nAfter Filtering") +
  geom_hline(yintercept = median(seurat.obj.org2$nCount_RNA), linetype = 'dashed')
vln3.org2 <- VlnPlot(seurat.obj.org2, features = c("percent.mt")) + labs(title = "Mitochondrial Content (%)\nAfter Removing Unwanted Cells") +
  geom_hline(yintercept = median(seurat.obj.org2@meta.data[["percent.mt"]]), linetype = 'dashed')

### Plot the final gene level UMAP----------------------------------
umap_title <- paste0("Org2 Sample\nResolution : ", cluster_resolution.org1, ", Dimensions: ", npc)
baseUMAPplot.org1 <- DimPlot(seurat.obj.org1, reduction = "umap", label = TRUE, ) + 
  labs(color = "cluster \n(from PCA)", 
       title = umap_title)
umap_title <- paste0("Org3 Sample\nResolution : ", cluster_resolution.org2, ", Dimensions: ", npc)
baseUMAPplot.org2 <- DimPlot(seurat.obj.org2, reduction = "umap", label = TRUE, ) + 
  labs(color = "cluster \n(from PCA)", 
       title = umap_title)

baseUMAPplot.org1 | baseUMAPplot.org2

### Create UMAPs showing the number of reads and unique genes per cell-----------------
## Create a UMAP plot showing the Reads per cell
UMIcountUMAP.org1 <- FeaturePlot(seurat.obj.org1, reduction = "umap", features = 'nCount_RNA') +
  labs(color = "UMI count",title = 'Number of Reads per Cell') + 
  theme(text = element_text(size = 10))
UMIcountUMAP.org2 <- FeaturePlot(seurat.obj.org2, reduction = "umap", features = 'nCount_RNA') +
  labs(color = "UMI count",title = 'Number of Reads per Cell') + 
  theme(text = element_text(size = 10))
UMIcountUMAP.org1 | UMIcountUMAP.org2

## Create a UMAP plot showing the number of genes expressed per cell
GeneCountUMAP.org1 <- FeaturePlot(seurat.obj.org1, reduction = "umap", features = 'nFeature_RNA')+
  labs(color = str_wrap("Gene count",15),title = 'Number of Unique Genes per Cell') + 
  theme(text = element_text(size = 10))
GeneCountUMAP.org2 <- FeaturePlot(seurat.obj.org2, reduction = "umap", features = 'nFeature_RNA')+
  labs(color = str_wrap("Gene count",15),title = 'Number of Unique Genes per Cell') + 
  theme(text = element_text(size = 10))
GeneCountUMAP.org1 | GeneCountUMAP.org2

## Create a UMAP plot showing mito content per cell
mitoPercUMAP.org1 <- FeaturePlot(seurat.obj.org1, reduction = "umap", features = 'percent.mt')+
  labs(color = str_wrap("mito %",15),title = 'Mitochondrial content per Cell') + 
  theme(text = element_text(size = 10))
mitoPercUMAP.org2 <- FeaturePlot(seurat.obj.org2, reduction = "umap", features = 'percent.mt')+
  labs(color = str_wrap("mito %",15),title = 'Mitochondrial content per Cell') + 
  theme(text = element_text(size = 10))
mitoPercUMAP.org1 | mitoPercUMAP.org2

### Export the processed seurat objects for use in further scripts
saveRDS(seurat.obj.org1, file = "org2-gene-QCed.rds")
saveRDS(seurat.obj.org2, file = "org3-gene-QCed.rds")

### perform the same processing steps for isoform level seurat objs---------

# Make sure that the barcodes that remain match those in the filtered gene level seurat object
# Such that the cells filtered out at the gene level are also removed in the transcript level object
seurat.obj.isoLvl.org1 <- subset(seurat.obj.isoLvl.org1, cells = row.names(seurat.obj.org1@meta.data))
seurat.obj.isoLvl.org2 <- subset(seurat.obj.isoLvl.org2, cells = row.names(seurat.obj.org2@meta.data))

seurat.obj.isoLvl.org1 <- NormalizeData(object = seurat.obj.isoLvl.org1)
seurat.obj.isoLvl.org1 <- FindVariableFeatures(seurat.obj.isoLvl.org1, 
                                                selection.method = 'vst',
                                                nfeatures = 2000)
seurat.obj.isoLvl.org2 <- NormalizeData(object = seurat.obj.isoLvl.org2)
seurat.obj.isoLvl.org2 <- FindVariableFeatures(seurat.obj.isoLvl.org2, 
                                               selection.method = 'vst',
                                               nfeatures = 2000)
  
all.features.org1 <- rownames(seurat.obj.isoLvl.org1)
seurat.obj.isoLvl.org1 <- ScaleData(seurat.obj.isoLvl.org1, features = all.features.org1)
seurat.obj.isoLvl.org1 <- RunPCA(seurat.obj.isoLvl.org1, features = VariableFeatures(object = seurat.obj.isoLvl.org1))

all.features.org2 <- rownames(seurat.obj.isoLvl.org2)
seurat.obj.isoLvl.org2 <- ScaleData(seurat.obj.isoLvl.org2, features = all.features.org2)
seurat.obj.isoLvl.org2 <- RunPCA(seurat.obj.isoLvl.org2, features = VariableFeatures(object = seurat.obj.isoLvl.org2))
  
#transfer cluster labels from gene object to isoform level object
seurat.obj.isoLvl.org1$seurat_clusters <- Idents(seurat.obj.org1)
Idents(seurat.obj.isoLvl.org1) <- seurat.obj.isoLvl.org1$seurat_clusters
seurat.obj.isoLvl.org2$seurat_clusters <- Idents(seurat.obj.org2)
Idents(seurat.obj.isoLvl.org2) <- seurat.obj.isoLvl.org2$seurat_clusters
  
seurat.obj.isoLvl.org1 <- RunUMAP(seurat.obj.isoLvl.org1, dims = 1:npc)
seurat.obj.isoLvl.org2 <- RunUMAP(seurat.obj.isoLvl.org2, dims = 1:npc)

#View(seurat.obj.isoLvl.org1@meta.data)
#View(seurat.obj.org1@meta.data)
#View(seurat.obj.isoLvl.org2@meta.data)
#View(seurat.obj.org2@meta.data)

## Generate a UMAP plot to see how cells cluster (using isoform counts) and match up with gene level cluster labels
baseUMAP.isoLvl.org1 <- DimPlot(seurat.obj.isoLvl.org1, reduction = "umap", label = TRUE, ) + 
    labs(color = "clusters \n(transferred from\ngene level counts)", 
         title = "Cerebral Organoid Sample 2\nIsoform Lvl UMAP")
baseUMAP.isoLvl.org2 <- DimPlot(seurat.obj.isoLvl.org2, reduction = "umap", label = TRUE, ) + 
  labs(color = "clusters \n(transferred from\ngene level counts)", 
       title = "Cerebral Organoid Sample 3\nIsoform Lvl UMAP")
baseUMAP.isoLvl.org1 | baseUMAP.isoLvl.org2
baseUMAPplot.org1 | baseUMAP.isoLvl.org1
baseUMAPplot.org2 | baseUMAP.isoLvl.org2

## Create a UMAP plot showing the number of unique isoforms expressed per cell
IsoformCountUMAP.org1 <- FeaturePlot(seurat.obj.isoLvl.org1, reduction = "umap", features = 'nFeature_RNA')+
    labs(color = str_wrap("Isoform count",15),title = 'Number of Unique Isoforms per Cell (Org2)') + 
    theme(text = element_text(size = 10))
IsoReadsUMAP.org1 <- FeaturePlot(seurat.obj.isoLvl.org1, reduction = "umap", features = 'nCount_RNA')+
  labs(color = str_wrap("UMI count",15),title = 'Number of Reads per Cell\n(Org2 Isoform Lvl)') + 
  theme(text = element_text(size = 10))
IsoReadsUMAP.org1 | IsoformCountUMAP.org1

IsoformCountUMAP.org2 <- FeaturePlot(seurat.obj.isoLvl.org2, reduction = "umap", features = 'nFeature_RNA')+
  labs(color = str_wrap("Isoform count",15),title = 'Number of Unique Isoforms per Cell (Org3)') + 
  theme(text = element_text(size = 10))
IsoReadsUMAP.org2 <- FeaturePlot(seurat.obj.isoLvl.org2, reduction = "umap", features = 'nCount_RNA')+
  labs(color = str_wrap("UMI count",15),title = 'Number of Reads per Cell\n(Org3 Isoform Lvl)') + 
  theme(text = element_text(size = 10))
IsoReadsUMAP.org2 | IsoformCountUMAP.org2

## Export graphs into a pdf
isolvl.layout <- grid.arrange(baseUMAP.isoLvl.org1, baseUMAPplot.org1, 
                              IsoformCountUMAP.org1, IsoReadsUMAP.org1,
                                nrow = 2, ncol = 2,
                                top=textGrob("QC Script - Isoform-Lvl Info\nCerebral Organoid Sample 2",
                                             gp = gpar(fontsize = 15)))
ggsave("org2-isoLvl-Summary-QCscript.pdf", isolvl.layout, width = 15, height = 15)

isolvl.layout <- grid.arrange(baseUMAP.isoLvl.org2, baseUMAPplot.org2, 
                              IsoformCountUMAP.org2, IsoReadsUMAP.org2,
                              nrow = 2, ncol = 2,
                              top=textGrob("QC Script - Isoform-Lvl Info\nCerebral Organoid Sample 3",
                                           gp = gpar(fontsize = 15)))
ggsave("org3-isoLvl-Summary-QCscript.pdf", isolvl.layout, width = 15, height = 15)


# Export the QCed transcript-level seurat objects
saveRDS(seurat.obj.isoLvl.org1, file = "org2-iso-QCed.rds")
saveRDS(seurat.obj.isoLvl.org2, file = "org3-iso-QCed.rds")

### Put together a summary table of parameters and statistics----------------------------------------
filtered.summary.org1 <- rbind("Sample ID" = "6M_org2",
                          "Initial Cell Num" = initial.cell.num.org1,
                          "Cells after QC" = dim(seurat.obj.org1)[2],
                          "Min. Gene Threshold" = min.features.org1,
                          "Max. Gene Threshold" = max.features.org1,
                          "Min. Reads Threshold" = min.counts,
                          "Max. Reads Threshold" = max.counts.sample1,
                          "Mitochondrial Content Threshold (%)" = MT_threshold,
                          "Median Genes per Cell" = median(seurat.obj.org1$nFeature_RNA),
                          "Median Isoforms per Cell" = median(seurat.obj.isoLvl.org1$nFeature_RNA),
                          "Median Reads per Cell"= median(seurat.obj.org1$nCount_RNA),
                          "Median MT % before filter" = round(median.mito.content.before.org1, 2),
                          "Median MT % after filter" = round(median(seurat.obj.org1@meta.data[["percent.mt"]]), 2),
                          "NPCs (UMAP dimensions)" = npc,
                          "Cluster Resolution" = cluster_resolution.org1)

filtered.summary.org2 <- rbind("Sample ID" = "6M_org2",
                               "Initial Cell Num" = initial.cell.num.org2,
                               "Cells after QC" = dim(seurat.obj.org2)[2],
                               "Min. Gene Threshold" = min.features.org2,
                               "Max. Gene Threshold" = max.features.org2,
                               "Min. Reads Threshold" = min.counts,
                               "Max. Reads Threshold" = max.counts.sample2,
                               "Mitochondrial Content Threshold (%)" = MT_threshold,
                               "Median Genes per Cell" = median(seurat.obj.org2$nFeature_RNA),
                               "Median Isoforms per Cell" = median(seurat.obj.isoLvl.org2$nFeature_RNA),
                               "Median Reads per Cell"= median(seurat.obj.org2$nCount_RNA),
                               "Median MT % before filter" = round(median.mito.content.before.org2, 2),
                               "Median MT % after filter" = round(median(seurat.obj.org2@meta.data[["percent.mt"]]), 2),
                               "NPCs (UMAP dimensions)" = npc,
                               "Cluster Resolution" = cluster_resolution.org2)

### Compile and arrange all the main QC figures that were generated ---------------------------------------
final.layout.alt.org1 <- grid.arrange(tableGrob(filtered.summary.org1), baseUMAPplot.org1,
                                 elbow.plt.org1,
                                 vln1.org1, vln2.org1,
                                 vln3.org1, association.plt.after.org1,
                                 doublets.umap.org1, tableGrob(statsDoublets.org1),
                                 UMIcountUMAP.org1, GeneCountUMAP.org1,
                                 mitoPercUMAP.org1,
                                 nrow = 6, ncol = 2,
                                 top = textGrob("SCRIPT 2 (QC) SUMMARY (6M_Org2)\n",
                                                gp = gpar(fontsize = 20)))
ggsave("org2-QC-summary-main.pdf", final.layout.alt.org1, width = 20, height = 30)

final.layout.alt.org2 <- grid.arrange(tableGrob(filtered.summary.org2), baseUMAPplot.org2,
                                      elbow.plt.org2,
                                      vln1.org2, vln2.org2,
                                      vln3.org2, association.plt.after.org2,
                                      doublets.umap.org2, tableGrob(statsDoublets.org2),
                                      UMIcountUMAP.org2, GeneCountUMAP.org2,
                                      mitoPercUMAP.org2,
                                      nrow = 6, ncol = 2,
                                      top = textGrob("SCRIPT 2 (QC) SUMMARY (6M_Org3)\n",
                                                     gp = gpar(fontsize = 20)))
ggsave("org3-QC-summary-main.pdf", final.layout.alt.org2, width = 20, height = 30)



### Save all the cluster resolution figures into a PDF---------------------------------------
cluster.res.layout.org1 <- grid.arrange(cluster_resolution.figs.org1[[1]], cluster_resolution.figs.org1[[2]],
                                   cluster_resolution.figs.org1[[3]], cluster_resolution.figs.org1[[4]], 
                                   cluster_resolution.figs.org1[[5]], cluster_resolution.figs.org1[[6]],
                                   nrow = 2, ncol = 3,
                                   top=textGrob(paste0("QC Script Org2 Cluster Resolution Figs"),
                                                gp = gpar(fontsize = 15)))
ggsave("org2-cluster-res-figures.pdf", cluster.res.layout.org1, width = 20, height = 18)

cluster.res.layout.org2 <- grid.arrange(cluster_resolution.figs.org2[[1]], cluster_resolution.figs.org2[[2]],
                                        cluster_resolution.figs.org2[[3]], cluster_resolution.figs.org2[[4]], 
                                        cluster_resolution.figs.org2[[5]], cluster_resolution.figs.org2[[6]],
                                        nrow = 2, ncol = 3,
                                        top=textGrob(paste0("QC Script Org3 Cluster Resolution Figs"),
                                                     gp = gpar(fontsize = 15)))
ggsave("org3-cluster-res-figures.pdf", cluster.res.layout.org2, width = 20, height = 18)