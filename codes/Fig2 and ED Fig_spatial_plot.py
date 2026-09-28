import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import squidpy as sq
from scipy.stats import mannwhitneyu
import tangram as tg

sc.logging.print_versions()
sc.set_figure_params(facecolor='white', figsize=(8, 8))
sc.settings.verbosity = 3

adata= sc.read_h5ad('./Visium_dataset.h5ad')
adata.shape


sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)
adata = adata[adata.obs['pct_counts_mt'] < 20].copy()
print(f'#cells after MT filter: {adata.n_obs}')
sc.pp.filter_genes(adata, min_cells=10)

sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=2000)

sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

sq.pl.spatial_scatter(
    adata,
    library_id='spatial',
    shape=None,
    color=['Cell_type']
)

sc.pl.spatial(
    adata,
    color="TLS_score",
    library_id='spatial',
    img_key='hires',
    alpha_img=0.15,
    color_map='coolwarm',
    edges=False
)

features= ['Cytotoxicity_signature', 'Exhaustion_signature']

for f in features:

    p = mannwhitneyu(
        adata.obs.loc[adata.obs['Neu16_status']=='Proximal to Neu16', f],
        adata.obs.loc[adata.obs['Neu16_status']=='Distal to Neu16', f]
    ).pvalue

    sns.violinplot(
        data=adata.obs,
        x='Neu16_status',
        y=f,
        inner='box'
    )

    plt.title(f'{'***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'ns'}')
    plt.show()

tg.plot_cell_annotation_sc(
    adata,
    ['Cell_clusters'],
    perc=0.03,
)