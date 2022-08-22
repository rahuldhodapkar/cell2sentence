#!/usr/bin/env python
## load_data_tutorial.py
#
# Tutorial testing automated dataset integration
#

import cell2sentence as cs
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import umap
import igraph as ig
import numpy as np


import os
os.makedirs('./calc', exist_ok=True)


import sys
import os
from pathlib import Path
import re
import pandas as pd
from tqdm import tqdm

import anndata
import cell2sentence as cs
import scanpy as sc
import scanpy.external as sce

species_tags = ['human', 'mouse', 'chick', 'zebrafish']
data_dir = './data'
outpath = './calc'

def enum_csv_samples(data_dir, tag):
    files = Path("{}/{}".format(data_dir, tag)).glob('*')
    return [f for f in files if f.is_file() and not f.stem.startswith('.')]

samples_to_process = enum_csv_samples(data_dir, 'human')

def read_csv_sample(s):
    adata = sc.read_csv(s)
    return(adata.T)

read_func = read_csv_sample,

adata_objs = []
for s in tqdm(samples_to_process[0:2]):
    adata = read_csv_sample(s)
    adata_objs.append(adata)

adata_combined = anndata.concat(adata_objs, axis=0)
csdata_combined = cs.transforms.csdata_from_adata(adata_combined)


# perform standard anndata processing
adata = adata_combined


######
# (1) filter cells and genes
######
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-') 
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

csdata_filtered = cs.transforms.csdata_from_adata(adata)

######
# (2) total-count normalize and logarithmize the data
######

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata.raw = adata

# continue to follow tutorial
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')


# dim reduc and plotting
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)

sc.tl.umap(adata)
sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'])

adata.obs['sample'] = np.array([x[:2] for x in adata.obs_names], dtype=object)
sc.pl.umap(adata, color='sample')

###################################################################
## RUN Harmony Integration
###################################################################

adata_harmony = adata
sce.pp.harmony_integrate(adata_harmony, 'sample')

adata_harmony.obsm['X_pca'] = adata_harmony.obsm['X_pca_harmony']
sc.pp.neighbors(adata_harmony, n_neighbors=15, n_pcs=40)
sc.tl.umap(adata_harmony)
sc.pl.umap(adata_harmony, color='sample')

###################################################################
## RUN CSDATA WORKFLOW
###################################################################

csdata = csdata_filtered

dist = csdata.create_distance_matrix(dist_type='levenshtein', prefix_len=20)
csdata.create_knn_graph(k=15)

reducer = umap.UMAP(metric='precomputed', n_components=2)
embedding = reducer.fit_transform(dist)

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=[0 if "H9" in x else 1 for x in csdata.cell_names],
    s=1,
    alpha=0.3)

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=csdata.get_rank_data_for_feature('ITGAX', invert=True),
    cmap = cmap,
    plotnonfinite = True,
    s=1)
plt.show()

clustering = csdata.knn_graph.community_leiden(
    objective_function='modularity')

fig, ax = plt.subplots()
ax.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=clustering.membership,
    label=clustering.membership,
    plotnonfinite = True,
    s=1)
plt.show()

markers_1 = csdata.find_differential_features(ident_1 = 
    [i for i, x in enumerate(clustering.membership) if x == 1])

markers_1.sort_values(by=['p_val'], ascending=True)

###################################################################
## Compute Diversity Index
###################################################################

g_unadj = ig.Graph.Weighted_Adjacency(adata.obsp['connectivities']).as_undirected()
g_harmony = ig.Graph.Weighted_Adjacency(adata_harmony.obsp['connectivities']).as_undirected()
g_csdata = csdata.knn_graph

g_csdata.diversity([0 if "H9" in x else 1 for x in csdata.cell_names])



