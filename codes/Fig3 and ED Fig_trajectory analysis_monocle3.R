expression_matrix <- readRDS('./Neutrophil_matrix.rds')
cell_metadata <- readRDS('./Neutrophil_metadata.rds')
gene_annotation <- readRDS('./Neutrophil_clusters_annotation.rds')

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 50)

cds <- align_cds(cds, alignment_group = 'Cohort')

cds <- reduce_dimension(cds)

cds <- cluster_cells(cds)

cds <- learn_graph(cds)

plot_cells(cds,
           color_cells_by = 'Neutrophil_clusters',
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
           
cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = 'Pseudotime',
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)