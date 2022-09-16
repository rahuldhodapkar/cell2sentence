#!/usr/bin/env python
## load_data_tutorial.py
#
# Tutorial for loading data and generating cell sentences.
#

import cell2sentence as cs
import scanpy as sc
from gensim.models import Word2Vec
import gensim
import numpy as np
import scipy as sp
import umap

import os
os.makedirs('./calc', exist_ok=True)

wv = gensim.models.KeyedVectors.load('./calc/gene_names_human_marrow.wv')

# standard cell2sentence workflow:
adata = sc.datasets.pbmc3k()
csdata = cs.transforms.csdata_from_adata(adata)

# (3) generate sentence lists
sentences = csdata.create_sentence_lists()

# (5) generate cell embeddings, and compute umap from model
embeddings = np.array([wv.get_mean_vector(s) for s in sentences])

word2vec_dist_mat = sp.spatial.distance.cdist(embeddings, embeddings, metric='cosine')
reducer = umap.UMAP(metric='precomputed')
word2vec_umap = reducer.fit_transform(word2vec_dist_mat)

# (6) add back to anndata object as custom embedding
adata.obsm['word2vec'] = embeddings
adata.obsm['word2vec_umap'] = word2vec_umap

sc.pp.neighbors(adata, n_neighbors=10, use_rep='word2vec')
sc.tl.leiden(adata)

# normalize anndata object for visualization, not needed for clustering
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


# calculate the top genes
sc.tl.rank_genes_groups(adata, groupby='leiden')
sc.pl.rank_genes_groups(adata, groupby='leiden')

sc.pl.embedding(adata, 'word2vec_umap',
    color=['leiden', 'CST3', 'NKG7', 'PPBP'],
    legend_loc='on data')

#
marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']
sc.pl.heatmap(adata, var_names=marker_genes, groupby='leiden')

