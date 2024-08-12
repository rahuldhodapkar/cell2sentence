#!/usr/bin/env python
#
# Tutorial for loading data and generating cell sentences.
#

import cell2sentence as cs
import scanpy as sc
from gensim.models import Word2Vec
import gensim
import numpy as np
import scipy as sp

# (0) Get a standard anndata object ready for analysis with csdata
#   - perform QC, filter low quality cells / doublets

adata = sc.datasets.pbmc3k()
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

################################################################################
## Begin Tutorial
################################################################################
#
# This tutorial will walk you through the basics of analyzing the PBMC3k
# dataset used in Seurat and scanpy tutorials.  We will use a pretrained
# word2vec embedding to embed and visualize the data without additional
# processing steps required.
#

# (1) load the pretrained human bone marrow trained vectors
wv = gensim.models.KeyedVectors.load(
    './pretrained/gene_embeddings/bone_marrow/gene_names_human_marrow.wv')

# (2) create csdata object from adata
csdata = cs.transforms.csdata_from_adata(adata)

# (3) generate sentence lists
sentences = csdata.create_sentence_lists()

# (4) generate cell embeddings, and compute umap from model
embeddings = np.array([wv.get_mean_vector(s) for s in sentences])

# (5) add back to anndata object as custom embedding
adata.obsm['word2vec'] = embeddings

# (6) continue analysis using scanpy utilities
sc.pp.neighbors(adata, n_neighbors=20, use_rep='word2vec')
sc.tl.umap(adata)
sc.tl.leiden(adata)

# normalize anndata object for visualization, not needed for clustering
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# plot a few marker genes
sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'], legend_loc='on data')

# plot more marker genes
marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']
sc.pl.umap(adata, color=marker_genes, legend_loc='on data')
