"""
Utility functions for cell2sentence
"""

#
# @author Rahul Dhodapkar
#

from collections import OrderedDict
from scipy import sparse
from tqdm import tqdm
from sklearn.utils import shuffle
from datasets import Dataset, DatasetDict

import numpy as np


def example_function(x):
    """
    Example function demonstrating how to import utilities. Returns 1
    """
    return 1


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
        print(
            (
                "WARN: more variables ({}) than observations ({})... "
                + "did you mean to transpose the object (e.g. adata.T)?"
            ).format(len(adata.var), len(adata.obs)),
            file=sys.stderr,
        )

    vocabulary = OrderedDict()
    gene_sums = np.ravel(np.sum(adata.X > 0, axis=0))

    for i, name in enumerate(adata.var_names):
        vocabulary[name.upper()] = gene_sums[i]  # keys are all uppercase gene names

    return vocabulary


def generate_sentences(adata, vocab, delimiter=' ', random_state=42):
    """
    Transform expression matrix to sentences. Sentences contain gene "words"
    denoting genes with non-zero expression. Genes are ordered from highest
    expression to lowest expression.

    Arguments:
        adata: an AnnData object to generate cell sentences from. Expects that
            `obs` correspond to cells and `vars` correspond to genes.
        vocab: an OrderedDict which as feature names as keys and counts as values
        random_state: sets the numpy random state for splitting ties
    Return:
        a `numpy.ndarray` of sentences, split by delimiter.
    """
    np.random.seed(random_state)

    if len(adata.var) > len(adata.obs):
        print(
            (
                "WARN: more variables ({}) than observations ({}), "
                + "did you mean to transpose the object (e.g. adata.T)?"
            ).format(len(adata.var), len(adata.obs)),
            file=sys.stderr,
        )

    mat = sparse.csr_matrix(adata.X)
    enc_map = list(vocab.keys())

    sentences = []
    for i in tqdm(range(mat.shape[0])):
        # for row i, [indptr[i]:indptr[i+1]] returns the indices of elements to take from 
        #  data and indices corresponding to row i
        cols = mat.indices[mat.indptr[i] : mat.indptr[i + 1]]
        vals = mat.data[mat.indptr[i] : mat.indptr[i + 1]]
        cols, vals = shuffle(cols, vals)
        sentences.append(delimiter.join([enc_map[x] for x in cols[np.argsort(-vals, kind="stable")]]))

    return sentences

def to_arrow_dataset(
    output_path: str, 
    cell_names: list, 
    sentences: list, 
    data_split_indices_dict: dict, 
    adata, 
    label_col_names: list
):
    """
    Write data represented by CSData to an arrow dataset.
    
    Arguments:
        output_path: save path where dataset should be written.
        cell_names: list of strings representing (unique) cell identifiers
        sentences: list of strings representing cell sentences
        data_split_indices_dict: dictionary of indices for train, val, and test sets
        adata: anndata.AnnData object
        label_col_names: list of column names in .obs DataFrame to save into dataset
                                along with cell sentences
    """
    data_dict = {
        "cell_name": cell_names,
        "cell_sentence": sentences,
    }
    if label_col_names is not None:
        for label_col in label_col_names:
            data_dict[label_col] = adata.obs[label_col].tolist()

    full_ds = Dataset.from_dict(data_dict)

    # Train/val/test split
    train_ds = full_ds.select(data_split_indices_dict["train"])
    val_ds = full_ds.select(data_split_indices_dict["val"])
    test_ds = full_ds.select(data_split_indices_dict["test"])

    # Create the DatasetDict with train, validation, and test splits
    ds_dict = DatasetDict({
        'train': train_ds,
        'validation': val_ds,
        'test': test_ds
    })
    
    # Save to disk
    ds_dict.save_to_disk(output_path)

