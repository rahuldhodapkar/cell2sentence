"""
Main data wrapper class definition
"""

#
# @authors: Rahul Dhodapkar, Syed Rizvi
#

# Python built-in libraries
import os

# Third-party libraries
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.utils import shuffle
from sklearn.model_selection import train_test_split
from datasets import load_from_disk
from tqdm import tqdm

# Local imports
# TODO: later on, change import back to .utils import ...
from cell2sentence.utils import generate_vocabulary, generate_sentences, to_arrow_dataset


class CSData():
    """
    Wrapper class to abstract different types of input data that can be passed
    in cell2sentence based workflows.
    """

    def __init__(self, vocab, feature_names, data_path, dataset_backend='arrow'):
        """
        Core constructor, CSData class contains a data path and format,
        and may also handle some buffering and data selection options.
        """
        self.vocab = vocab  # Ordered Dictionary: {gene_name: num_expressed_cells}
        self.feature_names = feature_names  # list of gene names
        self.data_path = data_path  # path to data file in arrow format
        self.dataset_backend = dataset_backend  # support plaintext and arrow

    @classmethod
    def from_adata(self, 
        adata, 
        save_dir: str, 
        save_name: str,
        label_col_names: list = None, 
        random_state: int = 42, 
        dataset_backend: str = 'arrow',
        sentence_delimiter: str = ' '
    ):
        """
        Create new CSData object from an anndata object.

        Arguments:
            adata: anndata.AnnData object to convert into a cell sentence dataset
            save_dir: directory where cell sentence dataset will be saved to disk
            save_name: name of folder to create storing cell sentence dataset (will be created)
            label_col_names: list of column names in .obs to save into dataset along with cell sentences
            random_state: random seed to control randomness
            dataset_backend: backend implementation for cell sentence datase
            sentence_delimiter: separator for cell sentence strings (default: ' ')
        """
        assert dataset_backend in ['arrow'], "C2S currently only supports arrow backend."
        
        # Create save directory
        save_path = os.path.join(save_dir, save_name)
        if not os.path.exists(save_path):
            os.makedirs(save_path, exist_ok=True)
        
        # Warn if var_names contains ensembl IDs instead of gene names.
        first_gene_name = str(adata.var_names[0])
        if "ENS" in first_gene_name:
            print(
                """WARN: adata.var_names seems to contain ensembl IDs rather than gene/feature names. 
                It is highly recommended to use gene names in cell sentences."""
            )

        # Create vocabulary and cell sentences based on adata object
        vocabulary = generate_vocabulary(adata)
        feature_names = list(vocabulary.keys())
        sentences = generate_sentences(adata, vocabulary, delimiter=sentence_delimiter)
        cell_names = adata.obs_names.tolist()

        # Train/val/test split here
        cell_indices_list = list(range(len(sentences)))
        train_and_val_indices, test_indices = train_test_split(cell_indices_list, test_size=0.1)
        train_indices, val_indices = train_test_split(train_and_val_indices, test_size=0.11)

        train_indices.sort()
        val_indices.sort()
        test_indices.sort()
        data_split_indices_dict = { "train": train_indices, "val": val_indices, "test": test_indices }
        
        # Save to disk
        if dataset_backend == "arrow":
            to_arrow_dataset(
                output_path=save_path, 
                cell_names=cell_names, 
                sentences=sentences,
                data_split_indices_dict=data_split_indices_dict,
                adata=adata,
                label_col_names=label_col_names
            )
        else:
            raise NotImplementedError("to_plain_text() function not yet implemented, please use arrow save format.")

        return self(
            vocab=vocabulary,
            feature_names=feature_names,
            data_path=save_path,
            dataset_backend='arrow',
        )

    def get_sentence_strings(self):
        """
        Helper function
        """
        ds_dict = load_from_disk(self.data_path)
        return {
            "train": ds_dict["train"]["cell_sentence"],
            "validation": ds_dict["validation"]["cell_sentence"],
            "test": ds_dict["test"]["cell_sentence"]
        }

    def to_plain_text(self, output_path):
        """
        Write data represented by CSData to a tab-separated plain text file.
        Arguments:
            output_path: a string representing the path to which the output file
                         should be written.
        """
        raise NotImplementedError("to_plain_text() function not yet implemented, please use arrow save format.")
    
    def __str__(self):
        """
        Summarize CSData object as string for debugging and logging.
        """
        return f"CSData Object; Path={self.data_path}, Format={self.dataset_backend}"
