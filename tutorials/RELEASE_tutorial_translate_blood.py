#!/usr/bin/env python
## translate_blood.py
#
# Tutorial for cross species translation of blood data
#

import scanpy as sc
import scipy as sp
import cell2sentence as cs
import gtfparse
import anndata
import gensim
import numpy as np
import pandas as pd
import plotnine as pn
import ot
from tqdm import tqdm

human_adata = sc.datasets.ebi_expression_atlas('E-MTAB-9221')
mouse_adata = sc.datasets.ebi_expression_atlas('E-CURD-117')

# convert gene ids to gene names
hu_df = gtfparse.read_gtf('/Users/rahuldhodapkar/Downloads/gencode.v41.chr_patch_hapl_scaff.annotation.gtf')
hu_ensembl_id_prefix_len = len('ENSG') + 11
hu_id2gene = {hu_df['gene_id'][i][:hu_ensembl_id_prefix_len]: hu_df['gene_name'][i] for i in range(hu_df.shape[0])}
human_adata.var_names = [hu_id2gene[x] if x in hu_id2gene else x for x in human_adata.var_names]
human_adata.var_names_make_unique()

mu_df = gtfparse.read_gtf('/Users/rahuldhodapkar/Downloads/gencode.vM10.annotation.gtf')
mu_ensembl_id_prefix_len = len('ENSMUSG') + 11
mu_id2gene = {mu_df['gene_id'][i][:mu_ensembl_id_prefix_len]: mu_df['gene_name'][i] for i in range(mu_df.shape[0])}
mouse_adata.var_names = [mu_id2gene[x] if x in mu_id2gene else x for x in mouse_adata.var_names]
mouse_adata.var_names_make_unique()

mouse_adata = sc.datasets.paul15() # for new mouse try
mouse_adata.obs['Sample Characteristic[organism]'] = 'Mus musculus'

# load pretrained gene vectors for translation and embedding
human_wv = gensim.models.KeyedVectors.load('./pretrained/gene_embeddings/bone_marrow/gene_names_human_marrow.wv')
mouse_wv = gensim.models.KeyedVectors.load('./pretrained/gene_embeddings/translate/bone_marrow/gene_names_mouse_to_human_marrow.wv')

human_csdata = cs.transforms.csdata_from_adata(human_adata)
mouse_csdata = cs.transforms.csdata_from_adata(mouse_adata)

human_embeddings = np.array([human_wv.get_mean_vector(s) for s in human_csdata.create_sentence_lists()])
mouse_embeddings = np.array([mouse_wv.get_mean_vector(s) for s in mouse_csdata.create_sentence_lists()])

# merge two anndata objects
human_adata.obsm['word2vec'] = human_embeddings
mouse_adata.obsm['word2vec'] = mouse_embeddings

# create unified cluster variable
human_adata.obs['celltype'] = human_adata.obs['Factor Value[inferred cell type - ontology labels]']
mouse_adata.obs['celltype'] = mouse_adata.obs['paul15_clusters']

###############################################################################
## Unify and Visualize
###############################################################################

adata = anndata.concat([human_adata, mouse_adata], join='outer')

# perform integrated analysis naturally with scanpy tools
sc.pp.neighbors(adata, n_neighbors=20, use_rep='word2vec')
sc.tl.umap(adata)

# plot combined representation
sc.pl.umap(adata, color='celltype', legend_loc='on data')

###############################################################################
## Generate Similarity Heatmap
###############################################################################

mouse_celltypes = mouse_adata.obs['celltype'].unique().tolist()
human_celltypes = human_adata.obs['celltype'].unique().tolist()
human_celltypes = human_celltypes[1:]

emd_dist_matrix = np.zeros((len(mouse_celltypes), len(human_celltypes)))
combined_dist_matrix = sp.spatial.distance.cdist(
    adata.obsm['word2vec'], adata.obsm['word2vec'], metric='cosine')

for i, c_from in enumerate(tqdm(mouse_celltypes)):
    for j, c_to in enumerate(human_celltypes):
        from_idxs = np.array(adata.obs['celltype']==c_from)
        to_idxs = np.array(adata.obs['celltype']==c_to)
        M = combined_dist_matrix[from_idxs,:][:,to_idxs]
        a = ot.unif(M.shape[0])
        b = ot.unif(M.shape[1])
        emd_dist_matrix[i,j] = ot.emd2(a,b,M, numItermax=1e7)

# generate heatmap
heatmap_df = pd.DataFrame({
    'from': pd.Categorical(np.concatenate(
        [[x] * len(human_celltypes) for x in mouse_celltypes]),
                           categories=mouse_celltypes),
    'to': pd.Categorical(human_celltypes * len(mouse_celltypes),
                           categories=human_celltypes),
    'affinity': np.ravel(1/emd_dist_matrix)
})

(pn.ggplot(heatmap_df, pn.aes(y='from', x='to', fill='affinity')) +
    pn.geom_tile() +
    pn.scales.scale_fill_cmap(cmap_name='magma') +
    pn.theme_classic())
