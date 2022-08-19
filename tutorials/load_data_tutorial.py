#!/usr/bin/env python
## load_data_tutorial.py
#
# Tutorial for loading data and generating cell sentences.
#

import cell2sentence as cs
import scanpy as sc
import matplotlib.pyplot as plt
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
dist = csdata.distance_matrix(dist_type='jaro', prefix_len=25)

# (4) compute UMAP embedding from distance matrix for visualization
reducer = umap.UMAP(metric='precomputed', n_components=2)
embedding = reducer.fit_transform(dist)

plt.scatter(
    embedding[:, 0],
    embedding[:, 1])


