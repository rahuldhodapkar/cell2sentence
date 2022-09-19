#!/usr/bin/env python
#
# Generate Mapping for two word vector files
#

import cell2sentence as cs
import scanpy as sc
from gensim.models import Word2Vec
import gensim
import numpy as np
import scipy as sp
import umap
import ot
from tqdm import tqdm

import os
os.makedirs('./calc/translate/retina', exist_ok=True)

retina_species_pairs = [
    ('mouse', 'human'),
    ('macaque', 'human'),
    ('zebrafish', 'human'),
    ('chick', 'human')
]

for src, tgt in retina_species_pairs:
    print("processing {} -> {}".format(src, tgt))
    #
    wv_species1 = gensim.models.KeyedVectors.load('./pretrained/gene_embeddings/retina/gene_names_{}_retina.wv'.format(src))
    wv_species2 = gensim.models.KeyedVectors.load('./pretrained/gene_embeddings/retina/gene_names_{}_retina.wv'.format(tgt))
    #
    homologs = []
    homolog_ixs = []
    for i, x in enumerate(tqdm(wv_species1.index_to_key)):
        for j, y in enumerate(wv_species2.index_to_key):
            if x.upper() == y.upper():
                homologs.append( (x,y) )
                homolog_ixs.append( (i,j) )
    #
    emb1 = wv_species1.vectors
    emb2 = wv_species2.vectors
    #
    words1 = wv_species1.index_to_key
    words2 = wv_species2.index_to_key
    #
    emb1 -= emb1.mean(axis=0)
    emb2 -= emb2.mean(axis=0)
    # and scale
    emb1 /= np.linalg.norm(emb1, axis=1)[:,None]
    emb2 /= np.linalg.norm(emb2, axis=1)[:,None]
    #
    M = np.ones((len(words1), len(words2)))
    for i, j in homolog_ixs:
        M[i, j] = 0
    #
    C1 = sp.spatial.distance.cdist(emb1, emb1, metric='cosine')
    C2 = sp.spatial.distance.cdist(emb2, emb2, metric='cosine')
    #
    C1 /= C1.mean()
    C2 /= C2.mean()
    #
    p = ot.unif(emb1.shape[0])
    q = ot.unif(emb2.shape[0])
    #
    gw0, log0 = ot.gromov.fused_gromov_wasserstein(
        M, C1, C2, p, q, loss_fun='square_loss', alpha=0.5, verbose=True, log=True)
    #
    gw = gw0 / np.reshape(np.sum(gw0, axis=1), newshape=(gw0.shape[0], 1))
    #
    # create new vectors and mappings for species 1 in species 2
    mapped_vectors = np.eye(len(wv_species1.index_to_key), 
                            len(wv_species1.index_to_key)) @ gw @ wv_species2.vectors
    #
    mapped_wv = gensim.models.KeyedVectors(vector_size=200)
    mapped_wv.add_vectors(
        words1,
        mapped_vectors
    )
    mapped_wv.save('./calc/translate/retina/gene_names_{}_to_{}_retina.wv'.format(src, tgt))




print("All done!")
