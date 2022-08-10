#!/usr/bin/env python
#
# Transform data from common structures used in the analysis of single-cell and
# single-nucleus RNA sequencing data, to cell sentences that can be used as
# input for natural language processing tools.
#
# @author Rahul Dhodapkar <rahul.dhodapkar@yale.edu>
#

from tqdm import tqdm
import os
import numpy as np
from scipy import sparse
import sys
from sklearn import model_selection

def train_test_validation_split(sentences, 
                                train_pct=0.6,
                                test_pct=0.2,
                                val_pct=0.2,
                                random_state=42):
    """
    Create train, test, and validation splits of the data given the supplied
    percentages with a specified random state for reproducibility.

    Arguments:
        sentences: an numpy.ndarray of sentences to be split.
        train_pct: Default = 0.6. the percentage of samples to assign to the training set.
        test_pct: Default = 0.2. the percentage of samples to assign to the test set.
        val_pct: Default = 0.2. the percentage of samples to assign to the validation set.
    Return:
        (train_sentences, test_sentences, val_sentences) split from the
        originally supplied sentences array.
    """
    if (train_pct + test_pct + val_pct != 1):
        raise ValueError(
            'train_pct = {} + test_pct = {} + val_pct = {} do not sum to 1.'.format(
            train_pct, test_pct, val_pct
        ))

    s1 = test_pct
    s2 = val_pct / (1 - test_pct)

    X = range(len(sentences))
    X_train, X_test = model_selection.train_test_split(
        X, test_size=s1, random_state=random_state)

    X_train, X_val = model_selection.train_test_split(
        X_train, test_size=s2, random_state=random_state)

    return (sentences[X_train], sentences[X_test], sentences[X_val])

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
    if (len(adata.var) > len(adata.obs)):
        print("WARN: more variables ({}) than observations ({}), did you mean to transpose the object (e.g. adata.T)?".format(
                len(adata.var), len(adata.obs)), file=sys.stderr)

    vocabulary = {}
    gene_sums = np.sum(adata.X > 0, axis=0)

    for i in range(len(adata.var_names)):
        vocabulary[adata.var_names[i]] = gene_sums[i]

    return(vocabulary)


def generate_sentences(adata, delimiter=" "):
    """
    Transform expression matrix to sentences. Sentences contain gene "words"
    denoting genes with non-zero expression. Genes are ordered from highest
    expression to lowest expression.

    Arguments:
        adata: an AnnData object to generate cell sentences from. Expects that
               `obs` correspond to cells and `vars` correspond to genes.
        delimiter: default = ' '. A token delimter for the generated sentences.
    Return:
        a `numpy.ndarray` of sentences, split by delimiter.
    """
    if (len(adata.var) > len(adata.obs)):
        print("WARN: more variables ({}) than observations ({}), did you mean to transpose the object (e.g. adata.T)?".format(
                len(adata.var), len(adata.obs)), file=sys.stderr)

    mat = sparse.csr_matrix(adata.X)

    sentences = []
    for i in tqdm(range(mat.shape[0])):
        cols = mat.indices[mat.indptr[i]:mat.indptr[i+1]]
        vals = mat.data[mat.indptr[i]:mat.indptr[i+1]]

        words = adata.var_names[cols[np.argsort(-vals)]]
        sentences.append(delimiter.join(words))

    return(np.array(sentences, dtype=object))
