#!/usr/bin/env python
#
# Cross species translation tutorial
#

import scanpy as sc
import cell2sentence as cs
from gensim.models import Word2Vec
import gensim
import numpy as np
import scipy as sp
import umap
import gtfparse

import mygene

####################################################################
## Run Human Model
####################################################################

# 596,816 cells, immune cell census

adata_human = adata = sc.datasets.ebi_expression_atlas("E-HCAD-4")
csdata = cs.transforms.csdata_from_adata(adata)
sentences = csdata.create_sentence_lists()
model = Word2Vec(sentences=sentences, vector_size=200, window=8, min_count=1, workers=4)

df = gtfparse.read_gtf('/Users/rahuldhodapkar/Downloads/gencode.v41.chr_patch_hapl_scaff.annotation.gtf')

ensembl_id_prefix_len = len('ENSG') + 11
id2gene = {df['gene_id'][i][:ensembl_id_prefix_len]: df['gene_name'][i] for i in range(df.shape[0])}

model.wv.save('./calc/id_human_marrow.wv')
# train the same model, but with gene symbols rather than gene IDs

gene_name_keys = [id2gene[x] if x in id2gene else x for x in model.wv.index_to_key]

wv_gene_names = gensim.models.KeyedVectors(vector_size=200)
wv_gene_names.add_vectors(
    gene_name_keys,
    model.wv.vectors
)

wv_gene_names.save('./calc/gene_names_human_marrow.wv')

model_gene_names.wv.ad

old_model = Word2Vec.load('./calc/human_marrow_old.model')

####################################################################
## Run Mouse Model
####################################################################

#adata_mouse = adata = sc.datasets.ebi_expression_atlas("E-ENAD-15")
#
# Download GTF from: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27/
#
df = gtfparse.read_gtf('/Users/rahuldhodapkar/Downloads/gencode.vM10.annotation.gtf')
gene2id = {df['gene_name'][i]: df['gene_id'][i] for i in range(df.shape[0])}
id2gene = {df['gene_id'][i][:ensembl_id_prefix_len]: df['gene_name'][i] for i in range(df.shape[0])}

#
# Download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4505404
# 
adata_mouse = adata = sc.read_h5ad('./data/tabula-muris-senis-droplet-official-raw-obj.h5ad')
adata = adata[adata.obs['tissue'] == 'Marrow',adata.var_names.isin(df['gene_name'])]

# get ensembl prefix length target
ensembl_id_prefix_len = len('ENSMUSG') + 11
mapped_ids = [gene2id[x][:ensembl_id_prefix_len] for x in adata.var_names]

adata.var_names = mapped_ids

csdata = cs.transforms.csdata_from_adata(adata)
sentences = csdata.create_sentence_lists()
model = Word2Vec(sentences=sentences, vector_size=200, window=8, min_count=1, workers=4)

model.wv.save('./calc/id_mouse_marrow.wv')

# now save another model with the gene names
gene_name_keys = [id2gene[x] if x in id2gene else x for x in model.wv.index_to_key]

wv_gene_names = gensim.models.KeyedVectors(vector_size=200)
wv_gene_names.add_vectors(
    gene_name_keys,
    model.wv.vectors
)
wv_gene_names.save('./calc/gene_names_mouse_marrow.wv')
