#!/usr/bin/env python
## load_data_tutorial.py
#
# Tutorial for loading data and generating cell sentences.
#

import cell2sentence as cs
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import umap
import igraph as ig

import os
os.makedirs('./calc', exist_ok=True)

# standard cell2sentence workflow:
adata = sc.read_10x_mtx(
    'tutorials/tutorialdata/filtered_gene_bc_matrices/hg19/',
    var_names='gene_symbols',
    cache=True)

# extract highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3')
adata_hvg = adata[:,adata.var.highly_variable]

# (2) convert to sentence format
csdata = cs.transforms.csdata_from_adata(adata_hvg)

# (3) generate edit distance matrix
dist = csdata.create_distance_matrix(dist_type='damerau_levenshtein', prefix_len=20)

# (4) create weighted k nearest neighbors graph
csdata.create_knn_graph(k=15)
clustering = csdata.knn_graph.community_multilevel()

c1_diff = csdata.find_differential_features(
    ident_1 = [i for i, x in enumerate(clustering.membership) if x == 1]
)
c1_diff.sort_values(by=['p_val'], ascending=True)

# (5) compute UMAP embedding from distance matrix for visualization
reducer = umap.UMAP(metric='precomputed', n_components=2)
embedding = reducer.fit_transform(dist)

# (6) visualize clusters on embedding
scatterplot = plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=clustering.membership,
    plotnonfinite = True,
    s=1)
plt.legend(*scatterplot.legend_elements(),
            loc="lower right", title="Clusters")
plt.show()

# (7) identify differential genes for each cluster.
diff_df = csdata.find_differential_features(
    ident_1 = [i for i, x in enumerate(clustering.membership) if x == 1]
)


diff_df = csdata.find_differential_features(
    ident_1 = [i for i, x in enumerate(clustering.membership) if x == 1]
)

#clustering = csdata.knn_graph.community_leiden(
#    objective_function='modularity')
#clustering = csdata.knn_graph.community_spinglass()
#clustering = csdata.knn_graph.community_walktrap().as_clustering()


# (7) plot characteristic gene expression by rank
cmap = mpl.cm.get_cmap("viridis").copy()
cmap.set_bad('gray', 0.7)

scatterplot = plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=csdata.get_rank_data_for_feature('PPBP', invert=True),
    cmap = cmap,
    plotnonfinite = True,
    s=1)
plt.legend(*scatterplot.legend_elements(),
            loc="lower right", title="Clusters")
plt.show()

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=csdata.get_rank_data_for_feature('LYZ', invert=True),
    cmap = cmap,
    plotnonfinite = True,
    s=1.5)
plt.legend()
plt.show()

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=csdata.get_rank_data_for_feature('MS4A1', invert=True),
    cmap = cmap,
    plotnonfinite = True,
    s=1.5)
plt.legend(*scatterplot.legend_elements(),
            loc="lower right", title="Clusters")
plt.show()

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=csdata.get_rank_data_for_feature('FCER1A', invert=True),
    cmap = cmap,
    plotnonfinite = True,
    s=1.5)
plt.show()

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=csdata.get_rank_data_for_feature('CD8A', invert=True),
    cmap = cmap,
    plotnonfinite = True,
    s=1.5)
plt.show()

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=csdata.get_rank_data_for_feature('FCGR3A', invert=True),
    cmap = cmap,
    plotnonfinite = True,
    s=1.5)
plt.show()

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=csdata.get_rank_data_for_feature('CD8A', invert=True),
    cmap = cmap,
    plotnonfinite = True,
    s=1.5)
plt.show()


