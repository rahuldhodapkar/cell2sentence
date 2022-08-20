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

import os
os.makedirs('./calc', exist_ok=True)

# standard cell2sentence workflow:
adata = sc.read_10x_mtx(
    'tutorials/tutorialdata/filtered_gene_bc_matrices/hg19/',
    var_names='gene_symbols',
    cache=True)

# (2) convert to sentence format
csdata = cs.transforms.csdata_from_adata(adata)

# (3) generate edit distance matrix
dist = csdata.create_distance_matrix(dist_type='jaro', prefix_len=25)

# (4) compute UMAP embedding from distance matrix for visualization
reducer = umap.UMAP(metric='precomputed', n_components=2)
embedding = reducer.fit_transform(dist)

cmap = mpl.cm.get_cmap("jet").copy()
cmap.set_bad('gray', 0.7)

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=csdata.get_rank_data_for_feature('PPBP', invert=True),
    cmap = cmap,
    plotnonfinite = True,
    s=1)
plt.show()

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=csdata.get_rank_data_for_feature('CD3E', invert=True),
    cmap = cmap,
    plotnonfinite = True,
    s=1.5)
plt.show()

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=csdata.get_rank_data_for_feature('MS4A1', invert=True),
    cmap = cmap,
    plotnonfinite = True,
    s=1.5)
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




