import numpy as np
from matplotlib import pyplot as plt
from pyslingshot import Slingshot
import scanpy as sc

adata= sc.read_h5ad('./Neutrophil.h5ad')

start_node= 'Neu_01_CD34'

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
custom_xlim = (-12, 12)
custom_ylim = (-12, 12)
# plt.setp(axes, xlim=custom_xlim, ylim=custom_ylim)

slingshot = Slingshot(adata, celltype_key="Neu_clusters", obsm_key="X_umap", start_node=start_node, is_debugging="verbose")

slingshot.fit(num_epochs=1, debug_axes=axes)

fig, axes = plt.subplots(ncols=2, figsize=(12, 4))
axes[0].set_title("Clusters")
axes[1].set_title("Pseudotime")
slingshot.plotter.curves(axes[0], slingshot.curves)
slingshot.plotter.clusters(axes[0], labels=np.arange(slingshot.num_clusters), s=4, alpha=0.5)
slingshot.plotter.clusters(axes[1], color_mode="pseudotime", s=5)