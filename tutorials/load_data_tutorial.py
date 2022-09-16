#!/usr/bin/env python
## load_data_tutorial.py
#
# Tutorial for loading data and generating cell sentences.
#

import cell2sentence as cs
import scanpy as sc
from gensim.models import Word2Vec
import numpy as np
import scipy as sp
import umap

import os
os.makedirs('./calc', exist_ok=True)

# standard cell2sentence workflow:
adata = sc.read_10x_mtx(
    'tutorials/tutorialdata/filtered_gene_bc_matrices/hg19/',
    var_names='gene_symbols',
    cache=True)

# follow the standard workflow
adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3')

csdata = cs.transforms.csdata_from_adata(adata[:, adata.var.highly_variable])

# continue standard normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]


# (3) generate sentence lists
sentences = csdata.create_sentence_lists()

# (4) train word2vec model
model = Word2Vec(sentences=sentences, vector_size=400, window=5, min_count=1, workers=4)

# we can now inspect the gene embeddings obtained
# 
model.wv['CD8A']
model.wv.most_similar('CD8B', topn=10)

# (5) generate cell embeddings, and compute umap from model
embeddings = np.array(
    [np.mean([model.wv[x] for x in s], axis=0) for s in sentences]
)

word2vec_dist_mat = sp.spatial.distance.cdist(embeddings, embeddings, metric='cosine')
reducer = umap.UMAP(metric='precomputed')
word2vec_umap = reducer.fit_transform(word2vec_dist_mat)

# (6) add back to anndata object as custom embedding
adata.obsm['word2vec'] = embeddings
adata.obsm['word2vec_umap'] = word2vec_umap

# can analyze as 
marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']

sc.pl.embedding(adata, 'word2vec_umap', color=['CD9'])
sc.pl.embedding(adata, 'word2vec_umap', color=['CST3', 'NKG7', 'PPBP'])

