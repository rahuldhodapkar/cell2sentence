"""
Data transformations and pre-processing for interpreting cells as sentences.

Transform data from common structures used in the analysis of single-cell and
single-nucleus RNA sequencing data, to cell sentences that can be used as
input for natural language processing tools.
"""

#
# @author Rahul Dhodapkar <rahul.dhodapkar@yale.edu>
#

import sys
from collections import OrderedDict
from itertools import chain

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.utils import shuffle

from tqdm import tqdm

from .csdata import CSData

def generate_vocabulary(adata):
    """
    Create a vocabulary dictionary, where each key represents a single gene
    token and the value represents the number of non-zero cells in the provided
    count matrix.

    Arguments:
        adata: an AnnData object to generate cell sentences from. Expects that
               `obs` correspond to cells and `vars` correspond to genes.
    Return:
        a dictionary of gene vocabulary
    """
    if len(adata.var) > len(adata.obs):
        print(("WARN: more variables ({}) than observations ({})... " +
              "did you mean to transpose the object (e.g. adata.T)?").format(
            len(adata.var), len(adata.obs)), file=sys.stderr)

    vocabulary = OrderedDict()
    gene_sums = np.ravel(np.sum(adata.X > 0, axis=0))

    for i, name in enumerate(adata.var_names):
        vocabulary[name] = gene_sums[i]

    return vocabulary


def generate_sentences(adata, prefix_len=None, random_state=42):
    """
    Transform expression matrix to sentences. Sentences contain gene "words"
    denoting genes with non-zero expression. Genes are ordered from highest
    expression to lowest expression.

    Arguments:
        adata: an AnnData object to generate cell sentences from. Expects that
               `obs` correspond to cells and `vars` correspond to genes.
        random_state: sets the numpy random state for splitting ties
    Return:
        a `numpy.ndarray` of sentences, split by delimiter.
    """
    np.random.seed(random_state)

    if len(adata.var) > len(adata.obs):
        print(("WARN: more variables ({}) than observations ({}), " +
               "did you mean to transpose the object (e.g. adata.T)?").format(
            len(adata.var), len(adata.obs)), file=sys.stderr)

    mat = sparse.csr_matrix(adata.X)

    sentences = []
    for i in tqdm(range(mat.shape[0])):
        cols = mat.indices[mat.indptr[i]:mat.indptr[i + 1]]
        vals = mat.data[mat.indptr[i]:mat.indptr[i + 1]]

        cols, vals = shuffle(cols, vals)

        sentences.append(
            ''.join([chr(x) for x in cols[np.argsort(-vals, kind='stable')]]))

    if prefix_len is not None:
        sentences = [s[:prefix_len] for s in sentences]

    return np.array(sentences, dtype=object)


def csdata_from_adata(adata, prefix_len=None, random_state=42):
    """
    Generate a CSData object from an AnnData object.

    Arguments:
        adata: an AnnData object to generate cell sentences from. Expects that
               `obs` correspond to cells and `vars` correspond to genes.
        prefix_len: consider only rank substrings of length prefix_len
        random_state: sets the numpy random state for splitting ties
    Return:
        a CSData object containing a vocabulary, sentences, and associated name data.
    """
    return CSData(
        vocab=generate_vocabulary(adata),
        sentences=generate_sentences(adata,
            prefix_len=prefix_len, random_state=random_state),
        cell_names=adata.obs_names,
        feature_names=adata.var_names
    )


def merge_csdata(csdata_lst):
    """
    Merge two csdata objects, assumes that features with the same name are the
    same feature and will collapse them accordingly.

    Arguments:
        csdata_lst: list of csdata objects to merge
    Return:
        a merged CSData object.
    """
    merged_features = set(chain.from_iterable(
        [x.vocab.keys() for x in csdata_lst]))
    merged_vocab = OrderedDict()

    for f in merged_features:
        merged_vocab[f] = 0
        for csdata in csdata_lst:
            merged_vocab[f] += csdata.vocab[f] if f in csdata.vocab else 0

    feat_to_merged_enc = {k: i for i, k in enumerate(merged_vocab.keys())}

    enc_maps = []
    for csdata in csdata_lst:
        enc_maps.append(
            {i: feat_to_merged_enc[k]
                for i, k in enumerate(csdata.vocab.keys())}
        )

    merged_sentences = []
    for i, csdata in enumerate(csdata_lst):
        for sent in csdata.sentences:
            merged_sentences.append(
                ''.join([chr(enc_maps[i][ord(x)]) for x in sent])
            )

    merged_cell_names = pd.Index(
        chain.from_iterable([x.cell_names for x in csdata_lst]), dtype=object)
    merged_feature_names = pd.Index(merged_vocab.keys(), dtype=object)

    return CSData(
        vocab=merged_vocab,
        sentences=merged_sentences,
        cell_names=merged_cell_names,
        feature_names=merged_feature_names
    )
