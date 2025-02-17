# This script performs an integrated analysis of single-cell RNA sequencing data using the Seurat package.
# It includes the following steps:
# 1. Loading necessary libraries and setting options.
# 2. Defining paths for data input and output.
# 3. Reading CellRanger data for multiple samples and creating Seurat objects.
# 4. Performing data integration and normalization.
# 5. Running dimensionality reduction (PCA and UMAP) and clustering.
# 6. Visualizing the integrated data using various plots.
# 7. Inferring cell cycle stages and regressing out cell cycle effects.
# 8. Identifying variable features and marker genes for clusters.
# 9. Saving results and generating heatmaps for marker genes.
# 10. Assigning cell type identities to clusters and identifying differentially expressed genes across conditions.

library(Seurat)
library(cowplot)
library(dplyr)

data_dir <- "path-to-data"
output_dir <- "output-dir"
cc_genes_file <- "regev_lab_cell_cycle_genes.txt"

# Read CellRanger data
samples <- c("WT_R1", "WT_R2", "CRE_R1", "CRE_R2", "PRE_R1", "PRE_R2", "LEU_R1", "LEU_R2", "MLEU_R1", "MLEU_R2")
data_list <- lapply(samples, function(sample) {
  data <- Read10X(data.dir = file.path(data_dir, sample))
  colnames(data) <- paste(colnames(data), sample, sep = '_')
  return(data)
})
names(data_list) <- samples

# Create Seurat objects
seurat_objects <- lapply(names(data_list), function(sample) {
  seurat_obj <- CreateSeuratObject(data_list[[sample]], project = "lymphoma", min.cells = 3, min.features = 200)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  seurat_obj$batch <- substr(sample, 1, 3)
  seurat_obj$celloforigin <- sample
  seurat_obj$sequencing <- as.character(match(sample, samples))
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & percent.mt < 5)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  return(seurat_obj)
})
names(seurat_objects) <- samples

# Perform integration
immune.anchors <- FindIntegrationAnchors(object.list = seurat_objects, dims = 1:30)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:30)

# clean up
rm(list = samples)

# Perform an integrated analysis
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.4)

# Save RDS file
saveRDS(object = immune.combined, file = file.path(output_dir, 'seurat3_all_noTLEU.rds'))

# Visualization
DimPlot(immune.combined, reduction = "umap", group.by = "batch", pt.size = 0.1)
DimPlot(immune.combined, reduction = "umap", label = TRUE, pt.size = 0.1)
DimPlot(immune.combined, reduction = "umap", group.by = "celloforigin", pt.size = 0.1)
DimPlot(immune.combined, reduction = "umap", group.by = "sequencing", pt.size = 0.1)
FeaturePlot(immune.combined, features = c("percent.mt"), pt.size = 0.1)

# Cell cycle inference
cc.genes <- readLines(con = cc_genes_file)
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
immune.combined <- CellCycleScoring(object = immune.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
DimPlot(immune.combined, reduction = "umap", group.by = "Phase", pt.size = 0.1)

# Regress Cell Cycle
immune.combined$CC.Difference <- immune.combined$S.Score - immune.combined$G2M.Score
immune.combined <- ScaleData(immune.combined, vars.to.regress = c("CC.Difference"), features = rownames(immune.combined))
immune.combined <- RunPCA(immune.combined, features = VariableFeatures(immune.combined), nfeatures.print = 10)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.35)
DimPlot(immune.combined, reduction = "umap", label = TRUE, pt.size = 0.1)
FeaturePlot(immune.combined, features = c("Hmgb2", "Stmn1", "H2afz", "Top2a"), ncol = 2, pt.size = 0.1)

# Plot Variable features
immune.combined <- FindVariableFeatures(immune.combined, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(immune.combined), 10)
VariableFeaturePlot(immune.combined)
LabelPoints(plot = VariableFeaturePlot(immune.combined), points = top10, repel = TRUE)

# Save table with stats
cl <- as.character(unique(immune.combined[[]]$seurat_clusters))
ct <- as.character(unique(immune.combined[[]]$celloforigin))
M <- matrix(0, nrow = length(cl), ncol = length(ct))
for (i in 1:length(cl)) {
  for (j in 1:length(ct)) {
    M[i, j] <- length(which(immune.combined[[]]$seurat_clusters == cl[i] & immune.combined[[]]$celloforigin == ct[j]))
  }
}
M <- as.data.frame(M)
rownames(M) <- cl
colnames(M) <- ct
M <- cbind(M, rowSums(M))
M <- rbind(M, colSums(M))
M <- data.frame(cluster = c(cl, 'colSums'), M)
write.csv(M, file = file.path(output_dir, 'All_cell_stats.csv'), row.names = FALSE)

# Find markers genes for all clusters
DefaultAssay(immune.combined) <- "RNA"
Nc <- as.character(unique(immune.combined[[]]$seurat_clusters))
for (i in 1:length(Nc)) {
  cluster.markers <- FindMarkers(immune.combined, ident.1 = Nc[i], min.pct = 0.1, logfc.threshold = log10(1.2))
  write.csv(cluster.markers, file = file.path(output_dir, paste('allcells_MarkerGenesCluster', Nc[i], 'csv', sep = '.')))
}

# Heatmap of gene markers (20 per cluster)
pbmc.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.1, logfc.threshold = log10(1.2))
top20 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
DoHeatmap(immune.combined, features = top20$gene, assay = "integrated") + NoLegend()

# Assigning cell type identity to clusters
new.cluster.ids <- paste0('C', 0:(length(Nc) - 1))
names(new.cluster.ids) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, new.cluster.ids)
DimPlot(immune.combined, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()

# Identify differential expressed genes across conditions
immune.combined$celltype.condition <- paste(Idents(immune.combined), immune.combined$batch, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.condition"
N <- as.character(new.cluster.ids)
for (j in 1:length(N)) {
  comp1 <- rep(paste0(N[j], "_WT"), 4)
  comp2 <- paste0(N[j], c('_CRE', '_PRE', '_LEU', '_MLEU'))
  for (i in 1:length(comp1)) {
    temp <- FindMarkers(immune.combined, ident.1 = comp1[i], ident.2 = comp2[i], verbose = FALSE, min.pct = 0.1, logfc.threshold = log10(1.2), min.cells.group = 1)
    write.csv(temp, file = file.path(output_dir, paste('DiffExp', comp1[i], 'vs', comp2[i], 'csv', sep = '.')))
  }
}

