#!/usr/bin/env Rscript

peptides_table <- "xxxx.txt"
labels_to_plot <- c( "DNMT3A", "DNMT3B" )
setwd("/path_to_txt_file/")

my_packages <- c("ggplot2", "ggrepel")
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]
if(length(not_installed)) install.packages(not_installed)

acetylated_peptides <- read.table( file = peptides_table, header=TRUE )
acetylated_peptides$significant <- ""
acetylated_peptides$genelabels <- ""

# P < 0.05
acetylated_peptides$significant[ acetylated_peptides$Log10adjPVal < 1.30103 ] <- "FALSE"
acetylated_peptides$significant[ acetylated_peptides$Log10adjPVal > 1.30103 ] <- "TRUE"

for ( k in 1:length(labels_to_plot ) ){ acetylated_peptides$genelabels [ which(acetylated_peptides$GeneSymbol == labels_to_plot[k] ) ] <- labels_to_plot [k] }

names( acetylated_peptides )
head( acetylated_peptides )

library( ggplot2 )
library( ggrepel )
ggplot( data = acetylated_peptides, aes(x = log2FC, y = Log10adjPVal, colour = significant, label = genelabels )) +  
    geom_point() +
    geom_text_repel( aes( label = genelabels ),
                     min.segment.length = 0,
                     segment.color = 'grey50',
                     seed = 42,
                     box.padding = 0.5, 
                     force = 20,
                     point.padding = 0.5,
                     max.overlaps = Inf ) +
    theme_minimal() 

ggsave( 'volcano_aceltyl_peptides.pdf', height = 7, width = 8 )

