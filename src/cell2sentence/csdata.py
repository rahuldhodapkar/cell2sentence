"""
Main data wrapper class definition
"""

#
# @author Rahul Dhodapkar
#

import os

import numpy as np
import pandas as pd
from scipy import sparse
from tqdm import tqdm
from datasets import Dataset, DatasetDict

from .utils import example_function


class CSData():
    """
    Wrapper class to abstract different types of input data that can be passed
    in cell2sentence based workflows.
    """

    def __init__(self, vocab, feature_names, data_path, data_path_format='arrow'):
        """
        Core constructor, CSData class contains a data path and format,
        and may also handle some buffering and data selection options.
        """
        self.vocab = vocab  # Ordered Dictionary: {gene_name: num_expressed_cells}
        self.feature_names = feature_names  # list of gene names
        self.data_path = data_path  # path to data file in arrow format
        self.data_path_format = data_path_format  # support plaintext and arrow

    @classmethod
    def from_adata(cls, 
        adata, 
        save_path: str, 
        label_col_name: str = None, 
        random_state: int = 42, 
        data_path_format: str = 'arrow'
    ):
        """
        Create new CSData object from an anndata object
        """
        assert data_path_format in ['arrow'], "Unknown data path format."  # TODO: add support for plaintext
        # TODO: consider adding a check that var_names contains gene names, not ensembl IDs.

        vocabulary = self.generate_vocabulary(adata)
        feature_names = adata.var_names
        sentences = self.generate_sentences(adata)
        cell_names = adata.obs_names

        data_dict = {
            "cell_name": cell_names,
            "cell_sentence": sentences,
        }
        if label_col_name is not None:
            data_dict["label"] = adata.obs["label_col_name"].tolist()

        full_ds = Dataset.from_dict(data_dict)

        # Train/val/test split
        train_test_split = full_ds.train_test_split(test_size=0.1)
        train_val_split = train_test_split['train'].train_test_split(test_size=0.1111)  # 10% of 90%

        # Create the DatasetDict with train, validation, and test splits
        ds_dict = DatasetDict({
            'train': train_val_split['train'],
            'validation': train_val_split['test'],
            'test': train_test_split['test']
        })

        # Create save directory
        if not os.path.exists(save_path):
            os.makedirs(save_path, exist_ok=True)
        
        # Save to disk
        ds_dict.save_to_disk(save_path)

        return cls(
            vocab=vocabulary,
            feature_names=feature_names,
            data_path=save_path,
            data_path_format='arrow'
        )

    def to_plain_text(self, output_path):
        """
        Print data represented by CSData to a plain text file
        Arguments:
            output_path: a string representing the path to which the output file
                         should be written.
        """
        return None

    def __str__(self):
        """
        Summarize CSData object as string for debugging and logging.
        """
        return f"CSData Object; Path={self.data_path}, Format={self.data_path_format}"

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
            print(
                (
                    "WARN: more variables ({}) than observations ({}), "
                    + "did you mean to transpose the object (e.g. adata.T)?"
                ).format(len(adata.var), len(adata.obs)),
                file=sys.stderr,
            )

        mat = sparse.csr_matrix(adata.X)
        # Scipy sparse_csr: data contains nonzero elems, indices maps elem to column, and 
        #  indptr maps to row index (for row i, [indptr[i]:indptr[i+1]] returns the indices 
        #  of elements to take from data and indices corresponding to row i. Tells you what 
        #  to take from indices and data)
        sentences = []
        for i in tqdm(range(mat.shape[0])):  # For each cell
            cols = mat.indices[mat.indptr[i] : mat.indptr[i + 1]]
            vals = mat.data[mat.indptr[i] : mat.indptr[i + 1]]

            cols, vals = shuffle(cols, vals)

            sentences.append(
                "".join([chr(x) for x in cols[np.argsort(-vals, kind="stable")]])
            )

        if prefix_len is not None:
            sentences = [s[:prefix_len] for s in sentences]

        return np.array(sentences, dtype=object)
