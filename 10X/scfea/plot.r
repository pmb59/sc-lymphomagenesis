if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)

col_names <- readLines("balance_20221130-210135.csv", n = 1)
metabol <- strsplit(col_names, split=',')[[1]][2:71]


mat=read.csv("balance_20221130-210135.csv")
mat2 <- as.matrix(mat[,2:ncol(mat)])
rownames(mat2) <- mat$X
colnames(mat2) <- metabol

ha = rowAnnotation(cluster=as.character( read.csv("meta.csv")$x ) ,
                   condition=read.csv("batch.csv")$x,
                   col = list(
                       condition = c("WT" = "cyan", "CRE" = "blue", "PRE" = "yellow",
                        "LEU" = "red", "MLEU" = "magenta")) )

pdf('scfea_mouse_atlas_batch.pdf', height=14, width=12 )
ht <- Heatmap(matrix=mat2, name="balance",  # try scale(mat2)
              column_title = "metabolic states",
              #row_title = "cells",
              right_annotation = ha,
              show_row_names=FALSE,
              row_split=as.character( read.csv("batch.csv")$x) )

draw(ht)
dev.off()


pdf('scfea_mouse_atlas_batch_scaled.pdf', height=14, width=12 )
ht <- Heatmap(matrix=scale(mat2), name="scaled balance", 
              column_title = "metabolic states",
              #row_title = "cells",
              right_annotation = ha,
              show_row_names=FALSE,
              row_split=as.character( read.csv("batch.csv")$x) )

draw(ht)
dev.off()
