# ==============================================================================
# Script Name: Fi1_Construction_of_singlecell_atlas.py
# Description: Construction of single-cell atlas
# Author: Zhang Lab
# ==============================================================================

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.pyplot import rc_context
from scib_metrics.benchmark import BatchCorrection, Benchmarker, BioConservation

# Global publication-quality figure configuration
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "Arial"

# adata_filtered = sc.read_h5ad("path_to_data.h5ad")

# ==============================================================================
# Figure S1C. Single-cell integration benchmarking (scIB) framework
# ==============================================================================

bm = Benchmarker(
    adata_filtered,
    batch_key="Cohort",
    label_key="clusters",
    bio_conservation_metrics=BioConservation(),
    batch_correction_metrics=BatchCorrection(),
    embedding_obsm_keys=["Unintegrated", "Scanorama", "Harmony", "BBKNN"],
)
# bm.benchmark()


# ==============================================================================
# Figure S2A, B, and C. UMAP with Cluster Labels
# ==============================================================================

palette = {"cell_cluster": "#FFFFFF"}  # Replace with actual color mapping
cluster_order = ["cell_cluster"]  # Replace with actual cluster order
cluster2num = {c: i + 1 for i, c in enumerate(cluster_order)}

# Extract UMAP coordinates
umap_df = pd.DataFrame(
    adata_filtered.obsm["X_umap"],
    columns=["UMAP1", "UMAP2"],
    index=adata_filtered.obs_names,
)
umap_df["cluster"] = adata_filtered.obs["clusters"].values

# Calculate cluster centroids
centers = umap_df.groupby("cluster")[["UMAP1", "UMAP2"]].median()
coords = centers.loc[cluster_order].values.copy()

# Iterative repulsion algorithm to prevent label overlapping
min_dist, step_size, max_iter = 2.5, 0.02, 2000
for _ in range(max_iter):
    moved = False
    for i in range(len(coords)):
        for j in range(i + 1, len(coords)):
            dx, dy = coords[i, 0] - coords[j, 0], coords[i, 1] - coords[j, 1]
            dist = np.sqrt(dx**2 + dy**2)
            if dist < min_dist:
                direction = (
                    np.array([1.0, 0.0])
                    if dist == 0
                    else np.array([dx, dy]) / dist
                )
                coords[i] += direction * step_size
                coords[j] -= direction * step_size
                moved = True
    if not moved:
        break

centers_offset = centers.loc[cluster_order].copy()
centers_offset["UMAP1"], centers_offset["UMAP2"] = coords[:, 0], coords[:, 1]

# Plot UMAP
with rc_context({"figure.figsize": (5, 5)}):
    ax = sc.pl.umap(
        adata_filtered,
        color="clusters",
        show=False,
        size=1,
        palette=palette,
        legend_loc=None,
    )

# Draw cluster labels with circle backgrounds
for cluster, row in centers_offset.iterrows():
    num = cluster2num[cluster]
    ax.scatter(
        row["UMAP1"],
        row["UMAP2"],
        s=580,
        facecolor=palette[cluster],
        edgecolor="none",
        zorder=10,
    )
    ax.text(
        row["UMAP1"],
        row["UMAP2"],
        str(num),
        ha="center",
        va="center",
        fontsize=18,
        fontweight="bold",
        color="black",
        zorder=10,
    )

# Axis cleanup and export
ax.set_xlabel("")
ax.set_ylabel("")
ax.set_xticks([])
ax.set_yticks([])
ax.set_title("")
plt.savefig("clusters_UMAP.pdf", format="pdf", dpi=300, bbox_inches="tight")
plt.close()


# ==============================================================================
# Figure S2D. Dot plot with Marginal Cluster Labels
# ==============================================================================

marker_genes = ["TXK", "TCF7", "LEF1"]  # Replace with actual marker genes
custom_cmap = LinearSegmentedColormap.from_list(
    "blue_to_red", ["#0d4179", "#549fc9", "#f8f7f3", "#d05747", "#890b26"]
)

dp = sc.pl.dotplot(
    adata_filtered,
    marker_genes,
    groupby="clusters",
    categories_order=cluster_order,
    standard_scale="var",
    color_map=custom_cmap,
    figsize=(13.3, 3.5),
    show=False,
)

ax = dp["mainplot_ax"]
plt.gcf().canvas.draw()
renderer = plt.gcf().canvas.get_renderer()

# Draw precise marginal indicators next to Y-axis text
for label in ax.get_yticklabels():
    cluster = label.get_text()
    if cluster not in cluster2num:
        continue

    num = cluster2num[cluster]
    bbox = label.get_window_extent(renderer=renderer)

    # Calculate pixel-based position with a slight left shift
    x_disp = bbox.x0 - 12
    y_disp = (bbox.y0 + bbox.y1) / 2
    x_data, y_data = ax.transData.inverted().transform((x_disp, y_disp))

    # Background circle
    ax.scatter(
        x_data,
        y_data,
        s=170,
        color=palette[cluster],
        edgecolor="none",
        clip_on=False,
        zorder=9,
    )

    # Label text
    ax.text(
        x_data,
        y_data,
        str(num),
        ha="center",
        va="center",
        fontsize=9,
        fontweight="bold",
        color="black",
        zorder=11,
    )

# Italicize gene names on X-axis
for label in ax.get_xticklabels():
    label.set_fontstyle("italic")

plt.savefig("genes_dotplot.pdf", format="pdf", dpi=300, bbox_inches="tight")
# plt.show()
plt.close()