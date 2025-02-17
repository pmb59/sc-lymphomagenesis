
# This script processes RNA sequencing data using the SingleR package.
# It takes a single argument which specifies the cluster to be processed.
#
# Example usage: Rscript --vanilla SingleR.r 0

library(Seurat)
library(cowplot)
library(dplyr)
library(SummarizedExperiment)
library(ExperimentHub)
library(AnnotationHub)
library(SingleR)
library(ggplot2)
x <- c(
  "SingleCellExperiment", "scater", "scRNAseq", "scran", 
  "BiocGenerics", "ggplot2", "pheatmap", "grDevices"
)
y <- c(
  "Rcpp", "beachmat", "methods", "Matrix", "S4Vectors", 
  "DelayedArray", "DelayedMatrixStats", "BiocNeighbors", 
  "BiocParallel", "stats", "utils", "Rcpp", "ExperimentHub", 
  "testthat", "knitr", "rmarkdown", "BiocStyle", "beachmat"
)
lapply(x, require, character.only = TRUE)
lapply(y, require, character.only = TRUE)


args <- commandArgs(trailingOnly=TRUE)
CLUSTER <-  args[1] 

setwd("path-to-working-directory")


# Load references from Immgen and RNAseq data
load(file="immgen.rda")
load("RNAseqForSingleR2.Rdata")
# remove duplicated gene names
expr_references <- expr_references[ which( duplicated(expr_references$external_gene_name) == FALSE) ,]
rownames(expr_references) <- expr_references$external_gene_name
expr_references <- expr_references[,1:18]
expr_references <- log10(expr_references+1)   # log FPKMs as recommended by SingleR
head(expr_references)
temp2 <- as.data.frame(expr_references)
temp_row_names <-  rownames(temp2)
immgen$data <- as.matrix(temp2)
rownames(immgen$data) <- temp_row_names
rm( temp2,temp_row_names)

# read counts from one Seurat cluster
Bcells.combined <- readRDS(file ='seurat3_all_noTLEU.rds')
head( Bcells.combined[[]] )

WT_Bcells  <- SubsetData(object=Bcells.combined, subset.name = "seurat_clusters", accept.value = CLUSTER )  # get cells for one cluster

# new variable for plotting barplot
LabelsBarplot <-  unique(WT_Bcells[[]]$batch)

test_data <- as.matrix(GetAssayData(object = WT_Bcells, slot = "data"))
test_condition <-  WT_Bcells[[]]$batch
seurat_clusters  <-  WT_Bcells[[]]$seurat_clusters
rm(WT_Bcells,Bcells.combined,x,y)


refAnnot <- read.csv('annotation_RNAseq2.csv', header=TRUE)
findS <- c()
for ( m in 1:nrow(refAnnot)  ){
  findS[m] <-  which( colnames(immgen$data) == refAnnot$sample[m] ) 
}

df  <- refAnnot 
immgen$data <- immgen$data[, findS]
rm(findS,refAnnot)

# run SingleR
pred.hpca <- SingleR(test = test_data, ref = immgen$data, labels = df$celltype , fine.tune=TRUE)  # fine.tune important parameter
pred.hpca
table(pred.hpca$labels)
sort(table(pred.hpca$labels), decreasing=TRUE)

# Save scores
SCORES <- data.frame( scores=pred.hpca$scores , first.labels=pred.hpca$first.labels, labels=pred.hpca$labels , pruned.labels = pred.hpca$pruned.labels ) 
rownames(SCORES) <- rownames(pred.hpca)
head(SCORES )
write.csv( SCORES, paste0("Cluster_",CLUSTER, '_scores_labels_all_RNAseq.csv' ) ) 
rm (SCORES)

# Plots
dfbarplot <- data.frame(annotation=pred.hpca$labels, seurat_clusters=as.character(seurat_clusters), condition=test_condition )
dfbarplot$condition <- factor(dfbarplot$condition, levels = LabelsBarplot )

p <- ggplot(dfbarplot, aes(x = annotation, fill = condition)) + 
  geom_bar(color = "black") +
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle = -90, hjust = 0),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(
    title = paste("Cluster", CLUSTER, "Annotation Barplot"),
    fill = "Condition"
  )
ggsave( paste0("Cluster_",CLUSTER, '_barplot_all_RNAseq.pdf' ), height=6, width=6)

my_heatmap <- plotScoreHeatmap(
  pred.hpca, 
  show.labels = TRUE, 
  show.pruned = FALSE,  
  normalize = TRUE, 
  max.labels = 100,
  annotation_col = data.frame(
    condition = test_condition,
    seurat_clusters = as.character(seurat_clusters),
    row.names = rownames(pred.hpca)
  )
)
save_pheatmap_png(my_heatmap, "my_heatmap.png")

save_pheatmap_pdf <- function(x, filename, width=9, height=14) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(my_heatmap, paste0("Cluster_", CLUSTER, "_my_heatmap_all_RNAseq.pdf")  )

# unnormarlised
my_heatmap2 <- plotScoreHeatmap(
  pred.hpca, 
  show.labels = TRUE, 
  show.pruned = FALSE,  
  normalize = FALSE, 
  max.labels = 100,
  annotation_col = data.frame(
    condition = test_condition,
    seurat_clusters = as.character(seurat_clusters),
    row.names = rownames(pred.hpca)
  )
)
save_pheatmap_pdf(my_heatmap2, paste0("Cluster_", CLUSTER, "_my_heatmap_unnormalisedScores_all_RNAseq.pdf")  )

