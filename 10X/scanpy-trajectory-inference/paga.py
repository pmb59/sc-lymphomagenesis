import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
import bbknn

sc.settings.set_figure_params(dpi=250, frameon=False)
sc.settings.verbosity = 3

clusterid = '8'

print("Scanpy and PAGA analysis")

adata = sc.read("counts.csv",delimiter='\t', cache=False)

sc.pp.recipe_zheng17(adata)

sc.tl.pca(adata, svd_solver='auto', zero_center='None')

condition = pd.read_csv("condition.csv") 
adata.obs['condition'] =  condition.values
adata.obs['condition'] =  adata.obs['condition'].astype('category')

anno = pd.read_csv("meta.csv") 
adata.obs['cell_groups'] =  anno.values
adata.obs['cell_groups'] =  adata.obs['cell_groups'].astype('category')

replicate = pd.read_csv("replicate.csv") 
adata.obs['replicate'] =  replicate.values
adata.obs['replicate'] =  adata.obs['replicate'].astype('category')

sequencing = pd.read_csv("sequencing.csv") 
adata.obs['sequencing'] =  sequencing.values
adata.obs['sequencing'] =  adata.obs['sequencing'].astype('category')

# Batch correction
corrData = bbknn.bbknn(adata, batch_key ='replicate')

# drawing single-cell graph using layout 'fa'
sc.tl.draw_graph(adata, layout='fa')

# ForceAtlas2
sc.pl.draw_graph(adata, color='condition', legend_fontsize=7, legend_loc='on data', size=3, title='Force Atlas', frameon='True', save=f"_ForceAtlasplot_cluster{clusterid}.pdf")

# Optional: Denoising the graph -- Diffussion map
# this is done representing the data in the Diffusion space (and not in PCA space)
print("denoising..")
sc.tl.diffmap(adata, n_comps=15)
sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_diffmap', knn='True')

sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata, layout='fa',color='condition', legend_loc='on data', size=10, legend_fontsize=10,title='Denoised Force Atlas', frameon='True', save=f"_DenoisedForceAtlas_cluster{clusterid}.pdf")

print("Running PAGA..")
sc.tl.paga(adata, groups='condition' )
sc.pl.paga(adata, threshold=0.20, fontsize=9,color=['condition'], layout='fa', title='PAGA graph (t=0.20)',  frameon='True', save=f"_PAGA_T=0.20_cluster{clusterid}.pdf")

# Recomputing the embedding using PAGA-initialization (trajectory on PAGA embeddings instead of ForceAtlas )
sc.tl.draw_graph(adata, init_pos='paga' )
sc.pl.draw_graph(adata, color=['condition'], layout='fa', size=10, legend_loc='on data', projection='2d', legend_fontsize=9,edges=False,  frameon='True', title='FA on PAGA embedding', save=f"_recomputing_cluster{clusterid}.pdf")

sc.tl.paga(adata, groups='condition' )
sc.pl.paga_compare( adata, threshold=0.05, title=f"Cluster {clusterid}", legend_loc='right margin', right_margin=0.3, size=10, legend_fontsize=7, fontsize=7, frameon=True, edges=True, save=f"_FA_PAGAembeddings_and_PAGA_graph_cluster{clusterid}.pdf")

# Diffusion Pseudotime
adata.uns['iroot'] = np.flatnonzero(adata.obs['condition']  == 'WT')[0]
sc.tl.dpt(adata)
sc.pl.draw_graph(adata, size=9,color=['condition', 'dpt_pseudotime'], legend_loc='right margin', legend_fontsize=7, save=f"_pseudotime_cluster{clusterid}.pdf")
