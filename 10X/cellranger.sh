#!/bin/bash

export PATH=/path_to/cellranger-5.0.0:$PATH

cellranger count --id=wt --transcriptome=/path_to/reference-genomes/refdata-cellranger-mm10-1.2.0 \
    --fastqs=./ --sample=SIGAH1 --localcores=8 \
    --localmem=64
