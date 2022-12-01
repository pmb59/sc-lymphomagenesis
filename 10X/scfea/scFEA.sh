#!/bin/bash

# metabolic flux using single-cell rna-seq data

conda activate scfea

cd scFEA
conda install --file requirements
conda install pytorch torchvision -c pytorch
pip install --user magic-impute

python src/scFEA.py --input_dir input --res_dir output \
  --test_file counts.csv --moduleGene_file module_gene_complete_mouse_m168.csv \
  --sc_imputation True
  
