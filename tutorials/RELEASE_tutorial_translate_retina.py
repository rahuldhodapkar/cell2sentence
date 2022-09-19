#!/usr/bin/env python
## translate_retina.py
#
# Tutorial for cross species translation of retina data
#

import scanpy as sc
import cell2sentence as cs
import gtfparse
import anndata
import gensim
import numpy as np

human_retina_adata = sc.datasets.ebi_expression_atlas('E-GEOD-137537')
mouse_retina_adata = sc.datasets.ebi_expression_atlas('E-MTAB-9061')

# convert gene ids to gene names
hu_df = gtfparse.read_gtf('/Users/rahuldhodapkar/Downloads/gencode.v41.chr_patch_hapl_scaff.annotation.gtf')
hu_ensembl_id_prefix_len = len('ENSG') + 11
hu_id2gene = {hu_df['gene_id'][i][:hu_ensembl_id_prefix_len]: hu_df['gene_name'][i] for i in range(hu_df.shape[0])}
human_retina_adata.var_names = [hu_id2gene[x] if x in hu_id2gene else x for x in human_retina_adata.var_names]
human_retina_adata.var_names_make_unique()

mu_df = gtfparse.read_gtf('/Users/rahuldhodapkar/Downloads/gencode.vM10.annotation.gtf')
mu_ensembl_id_prefix_len = len('ENSMUSG') + 11
mu_id2gene = {mu_df['gene_id'][i][:mu_ensembl_id_prefix_len]: mu_df['gene_name'][i] for i in range(mu_df.shape[0])}
mouse_retina_adata.var_names = [mu_id2gene[x] if x in mu_id2gene else x for x in mouse_retina_adata.var_names]
mouse_retina_adata.var_names_make_unique()

# load pretrained gene vectors for translation and embedding
human_wv = gensim.models.KeyedVectors.load('./pretrained/gene_embeddings/retina/gene_names_human_retina.wv')
mouse_wv = gensim.models.KeyedVectors.load('./pretrained/translate/retina/gene_names_mouse_to_human_retina.wv')


human_csdata = cs.transforms.csdata_from_adata(human_retina_adata)
mouse_csdata = cs.transforms.csdata_from_adata(mouse_retina_adata)

human_embeddings = np.array([human_wv.get_mean_vector(s) for s in human_csdata.create_sentence_lists()])
mouse_embeddings = np.array([mouse_wv.get_mean_vector(s) for s in mouse_csdata.create_sentence_lists()])

# merge two anndata objects
human_retina_adata.obsm['word2vec'] = human_embeddings
mouse_retina_adata.obsm['word2vec'] = mouse_embeddings

adata = anndata.concat([human_retina_adata, mouse_retina_adata], join='outer')

# perform integrated analysis naturally with scanpy tools
sc.pp.neighbors(adata, n_neighbors=20, use_rep='word2vec')
sc.tl.umap(adata)

# plot combined representation
sc.pl.umap(adata, color='Sample Characteristic[organism]', legend_loc='on data')

sc.pl.umap(adata, color=['RHO', 'Rho'])

sc.pl.umap(adata, color=['GLUL', 'Glul'])

sc.pl.umap(adata, color=['CDH5', 'Cdh5'])

sc.pl.umap(adata, color='Factor Value[inferred cell type - ontology labels]', legend_loc='on data')


