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
import plotnine as pn
import pandas as pd

import os
os.makedirs('./calc', exist_ok=True)

# standard cell2sentence workflow:
adata = sc.read_10x_mtx(
    'tutorials/tutorialdata/filtered_gene_bc_matrices/hg19/',
    var_names='gene_symbols',
    cache=True)

# filter adata object
adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

#################
# extract highly variable genes and save to csdata
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3')
adata_hvg = adata[:,adata.var.highly_variable]
csdata = cs.transforms.csdata_from_adata(adata_hvg)
#################

# normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata)

sc.pl.umap(adata, color=['leiden', 'CST3', 'NKG7'])
new_cluster_names = [
    'CD4 T', 
    'CD14 Monocytes',
    'B', 
    'CD8 T',
    'FCGR3A Monocytes',
    'NK',
    'Dendritic',
    'Megakaryocytes']
adata.rename_categories('leiden', new_cluster_names)

plot_df = pd.DataFrame({
    'umap1': adata.obsm['X_umap'][:,0],
    'umap2': adata.obsm['X_umap'][:,1],
    'scanpy_clusters': adata._obs['leiden']
})
p = (pn.ggplot(plot_df, pn.aes(x='umap1', y='umap2', color='scanpy_clusters')) +
    pn.geom_point(size=1) + pn.theme_classic())
p.save('./calc/scanpy_umap_labeled.png', width=8, height=8)



dists_to_compare = ['levenshtein', 'damerau_levenshtein', 'zlib_ncd', 'jaro', 'jaro_winkler']

for distance_method in dists_to_compare:
    # (3) generate edit distance matrix
    dist = csdata.create_distance_matrix(dist_type=distance_method, prefix_len=20)
    #
    # (4) create weighted k nearest neighbors graph
    csdata.create_knn_graph(k=15)
    clustering = csdata.knn_graph.community_multilevel()
    #
    # (5) compute UMAP embedding from distance matrix for visualization
    reducer = umap.UMAP(metric='precomputed', n_components=2)
    embedding = reducer.fit_transform(dist)
    #
    plot_df = pd.DataFrame({
        'umap1': embedding[:, 0],
        'umap2': embedding[:, 1],
        'scanpy_clusters': adata._obs['leiden']
    })
    #
    p = (pn.ggplot(plot_df, pn.aes(x='umap1', y='umap2', color='scanpy_clusters')) +
        pn.geom_point(size=1) + pn.theme_classic())
    p.save('./calc/c2s_{}_umap_labeled.png'.format(distance_method), width=8, height=8)


# (6) visualize clusters on embedding
scatterplot = plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=adata._obs['leiden'],
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



