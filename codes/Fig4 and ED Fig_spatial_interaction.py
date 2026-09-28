import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import sys
import laris as la

# Colours
from matplotlib import cm
from matplotlib import colors, colorbar
cmap_own = cm.get_cmap('magma_r', 256)
newcolors = cmap_own(np.linspace(0,0.75 , 256))
Greys = cm.get_cmap('Greys_r', 256)
#newcolors[:1, :] = Greys(np.linspace(0.8125, 0.8725, 1))
newcolors[:10, :] = Greys(np.linspace(0.8125, 0.8725, 10))
pos_cmap = colors.ListedColormap(newcolors)

adata=sc.read('./Visium_dataset.h5ad')
adata.shape

lr_df=pd.read_csv('./human_lr_database.csv',index_col=0)

rows_keep=np.logical_and( lr_df['ligand'].isin(adata.var_names), lr_df['receptor'].isin(adata.var_names) )

lr_df=lr_df.loc[rows_keep].copy()

lr_adata=la.tl.prepareLRInteraction(
    adata, 
    lr_df, 
    number_nearest_neighbors=20,
    use_rep_spatial='X_spatial'
)

LARIS_variable_interactions, res_LARIS=la.tl.runLARIS(
    lr_adata,
    adata,
    use_rep='X_spatial',
    n_nearest_neighbors=20,
    random_seed=27,
    n_repeats=5,
    mu = 0.40, # larger value indicates more consideration for specificity
    sigma=100,
    remove_lowly_expressed=False,
    expressed_pct=0.1,
    n_cells_expressed_threshold=100,
    n_top_lr=lr_adata.shape[1],

    by_celltype=True,
    ### Parameters for cell type-level inference when by_celltype set to True
    groupby='Neu_clusters', # label to group by
    use_rep_spatial='X_spatial', # spatial coordinates
    mu_celltype=100, # higher value puts more emphasis on cell type specificity
    expressed_pct_celltype=0.1, # expression percentage cut off per cell type
    remove_lowly_expressed_celltype=False,
    mask_threshold=1e-6, # mask for cosg, lower value is less restrictive
    n_neighbors_permutation=30,
    score_threshold= 1e-10,
    spatial_weight = 3.0,
)

x_width=adata.obsm['X_spatial'][:,0].max()-adata.obsm['X_spatial'][:,0].min()
y_width=adata.obsm['X_spatial'][:,1].max()-adata.obsm['X_spatial'][:,1].min()

plt.rcParams['figure.figsize'] = 2.5, 2.5*y_width/x_width

sc.pl.embedding(
    lr_adata,
    basis='X_spatial',
    color=['CCL19::CCR7','CD40LG::CD40','CXCL13::CXCR5'],
    cmap=pos_cmap,
    ncols=4,
    # size=120,
    frameon=False)
plt.rcParams['figure.figsize'] = 4, 4












