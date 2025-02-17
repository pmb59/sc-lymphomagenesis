# This script performs an integration analysis of single-cell RNA sequencing data from mouse and human samples.
# It reads the data, processes it, and integrates it using the Seurat package.
# The script includes steps for data normalization, variable feature selection, scaling, PCA, UMAP, and clustering.
# It also generates various visualizations to explore the integrated data.
# Finally, it saves the integrated Seurat object to an RDS file for future use.

library(Seurat)
library(cowplot)
library(dplyr)
library(biomaRt)

data_dir <- "path-to-data"
mouse_dir <- "path-to-mouse-data"
gene_file <- 'mouse_genes.tsv'
humanFL_file <- 'GSM3190076_humanFL5_matrixGeneCell.csv'

setwd( data_dir )

# Get Human Genes Equivalent
LPJ119_1.data <- Read10X(data.dir = file.path(data_dir, "LPJ119_1"), gene.column = 1)
HG <- LPJ119_1.data@Dimnames[[1]]

ensemblGRCh37 <- useMart(host = 'http://grch37.ensembl.org',
             biomart = 'ENSEMBL_MART_ENSEMBL',
             dataset = 'hsapiens_gene_ensembl')
geneinfo <- getBM(
    attributes = c(
        'ensembl_gene_id', 
        'mmusculus_homolog_orthology_type', 
        'mmusculus_homolog_associated_gene_name'
    ), 
    mart = ensemblGRCh37
)

new_dimnames <- HG
for (j in 1:length(HG)) {
  ix <- which(geneinfo$ensembl_gene_id == HG[j])
  if (length(ix) == 1) {
    if (geneinfo$mmusculus_homolog_orthology_type[ix] == "ortholog_one2one" & nchar(geneinfo$mmusculus_homolog_associated_gene_name[ix]) > 0) {
        new_dimnames[j] <- geneinfo$mmusculus_homolog_associated_gene_name[ix]
    }
  }
}
iz <- which(duplicated(new_dimnames) == TRUE)
new_dimnames[iz] <- paste0(new_dimnames[iz], '.dup')
rm(LPJ119_1.data, HG, geneinfo, ix, iz, j, ensemblGRCh37)

# Read CellRanger mouse data
samples <- c("WT_R1", "WT_R2", "CRE_R1", "CRE_R2", "PRE_R1", "PRE_R2", "LEU_R1", "LEU_R2", "MLEU_R1", "MLEU_R2")
data_list <- lapply(samples, function(sample) {
  data <- Read10X(data.dir = file.path(mouse_dir, sample))
  colnames(data) <- paste(colnames(data), sample, sep = '_')
  data
})

# Read human data
new_samples <- c("LPJ119_1", "LPJ119_2", "LPJ128A_1", "LPJ128A_2", "LPJ128A_3", "LPM011_1", "LPM011_2", "LPM018A_1", "LPM018A_2", "LPM020_2", "LPM021_1", "LPM021_2", "Tonsil04_2", "Tonsil04_3", "Tonsil92_1", "Tonsil92_2")
new_data_list <- lapply(new_samples, function(sample) {
  data <- Read10X(data.dir = file.path(data_dir, sample), gene.column = 1)
  data@Dimnames[[1]] <- new_dimnames
  colnames(data) <- paste(colnames(data), sample, sep = '_')
  data
})

# Combine all data
all_data <- c(data_list, new_data_list)

# Create Seurat objects
seurat_objects <- lapply(all_data, function(data) {
  seurat_obj <- CreateSeuratObject(data, project = "lymphoma", min.cells = 3, min.features = 200)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & percent.mt < 5)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj
})

# Set metadata for Seurat objects
metadata <- data.frame(
  batch = c(rep("WT", 2), rep("CRE", 2), rep("PRE", 2), rep("LEU", 2), rep("MLEU", 2), rep("LPJ119", 2), rep("LPJ128A", 3), rep("LPM011", 2), rep("LPM018A", 2), "LPM020", rep("LPM021", 2), rep("Tonsil04", 2), rep("Tonsil92", 2)),
  celloforigin = c(samples, new_samples),
  species = c(rep("mouse", 10), rep("human", 16)),
  study = c(rep("Horton et al", 10), rep("Andor et al", 16))
)

for (i in 1:length(seurat_objects)) {
  seurat_objects[[i]]$batch <- metadata$batch[i]
  seurat_objects[[i]]$celloforigin <- metadata$celloforigin[i]
  seurat_objects[[i]]$species <- metadata$species[i]
  seurat_objects[[i]]$study <- metadata$study[i]
}

# Perform integration
immune.anchors <- FindIntegrationAnchors(object.list = seurat_objects, dims = 1:30)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:30)

# Total Number of cells after initial filtering and preprocessing
print(nrow(immune.combined[[]]))

# Perform an integrated analysis
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.4)

# Save RDS file
saveRDS(object = immune.combined, file = 'seurat3_mouse_human.rds')

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "batch", pt.size = 0.1)
p1
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, pt.size = 0.1)
p2
p3 <- DimPlot(immune.combined, reduction = "umap", group.by = "celloforigin", pt.size = 0.1)
p3
p4 <- DimPlot(immune.combined, reduction = "umap", group.by = "species", pt.size = 0.1)
p4
p4b <- FeaturePlot(immune.combined, features = c("percent.mt"), pt.size = 0.1)
p4b
p5 <- DimPlot(immune.combined, reduction = "umap", group.by = "study", pt.size = 0.1)
p5

# TSNE
immune.combined <- RunTSNE(object = immune.combined, reduction.use = "cca.aligned", dims = 1:30)
DimPlot(immune.combined, reduction = "tsne", label = TRUE, pt.size = 0.1)

# Feature plots
FeaturePlot(immune.combined, features = c("Bin1", "Cd36", "Cdk1", "Mcm7"), ncol = 2, pt.size = 0.1, cols = c("black", "yellow", 'red'))

# UMAP by condition
DimPlot(immune.combined, reduction = "umap", split.by = "batch", ncol = 2)

# UMAP by replicate
DimPlot(immune.combined, reduction = "umap", split.by = "celloforigin", ncol = 2)

# UMAP by species
DimPlot(immune.combined, reduction = "umap", split.by = "species", ncol = 2)

# UMAP by study
DimPlot(immune.combined, reduction = "umap", split.by = "study", ncol = 2)

# Differential expression tests on the “unintegrated” data
DefaultAssay(immune.combined) <- "RNA"

# Save Seurat file
saveRDS(immune.combined, file = 'seurat3_mouse_human.rds')

# Feature plots
FeaturePlot(immune.combined, features = c("Bin1", "Cd36", "Cdk1", "Mcm7"), ncol = 2, pt.size = 0.1, cols = c("black", "yellow", 'red'))
FeaturePlot(immune.combined, features = c("Aicda"), pt.size = 0.1)
FeaturePlot(immune.combined, features = c("Rgs13"), pt.size = 0.1)
FeaturePlot(immune.combined, features = c("S1pr2"), pt.size = 0.1)
FeaturePlot(immune.combined, features = c("Marcksl1"), pt.size = 0.1)
FeaturePlot(immune.combined, features = c("Sdc1"), pt.size = 0.1)
FeaturePlot(immune.combined, features = c("Xbp1"), pt.size = 0.1)
FeaturePlot(immune.combined, features = c("Tppp3"), pt.size = 0.1)
FeaturePlot(immune.combined, features = c("Cd36"), pt.size = 0.1)

