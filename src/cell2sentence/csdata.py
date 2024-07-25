"""
Main data wrapper class definition
"""

#
# @authors: Rahul Dhodapkar, Syed Rizvi
#

import os

import numpy as np
import pandas as pd
from scipy import sparse
from tqdm import tqdm
from sklearn.utils import shuffle
from sklearn.model_selection import train_test_split
from datasets import load_from_disk

from cell2sentence.utils import generate_vocabulary, generate_sentences, to_arrow_dataset


class CSData():
    """
    Wrapper class to abstract different types of input data that can be passed
    in cell2sentence based workflows.
    """

    def __init__(self, vocab, feature_names, data_path, data_split_indices_dict, data_path_format='arrow'):
        """
        Core constructor, CSData class contains a data path and format,
        and may also handle some buffering and data selection options.
        """
        self.vocab = vocab  # Ordered Dictionary: {gene_name: num_expressed_cells}
        self.feature_names = feature_names  # list of gene names
        self.data_path = data_path  # path to data file in arrow format
        self.data_split_indices_dict = data_split_indices_dict  # train, val, and test split indices
        self.data_path_format = data_path_format  # support plaintext and arrow

    @classmethod
    def from_adata(self, 
        adata, 
        save_path: str, 
        label_col_name: str = None, 
        random_state: int = 42, 
        data_path_format: str = 'arrow',
        delimiter: str = ' '
    ):
        """
        Create new CSData object from an anndata object
        """
        assert data_path_format in ['arrow'], "Unknown data path format."  # TODO: add support for plaintext
        
        # Create save directory
        if not os.path.exists(save_path):
            os.makedirs(save_path, exist_ok=True)
        
        # Check that var_names contains gene names, not ensembl IDs.
        first_gene_name = str(adata.var_names[0])
        if "ENS" in first_gene_name:
            print("WARN: adata.var_names seems to contain ensembl IDs rather than gene/feature names. It is highly recommended to use gene names in cell sentences.")

        # Create vocabulary and cell sentences based on adata object
        vocabulary = generate_vocabulary(adata)
        feature_names = list(vocabulary.keys())
        sentences = generate_sentences(adata, vocabulary, delimiter=delimiter)
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
        if data_path_format == "arrow":
            to_arrow_dataset(
                output_path=save_path, 
                cell_names=cell_names, 
                sentences=sentences,
                data_split_indices_dict=data_split_indices_dict,
                label_col_name=label_col_name
            )
        else:
            raise NotImplementedError("to_plain_text() function not yet implemented, please use arrow save format.")

        return self(
            vocab=vocabulary,
            feature_names=feature_names,
            data_path=save_path,
            data_split_indices_dict=data_split_indices_dict,
            data_path_format='arrow',
        )

    def get_sentence_strings(self):
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
        return f"CSData Object; Path={self.data_path}, Format={self.data_path_format}"


# Debugging
if __name__ == "__main__":
    from pathlib import Path
    import scanpy as sc

    HERE = Path(__file__).parent
    adata = sc.read_csv(HERE / 'tests/small_data.csv').T
    csdata = CSData.from_adata(adata, save_path="/home/sr2464/palmer_scratch/C2S_Files_Syed/c2s_api_testing")
    cell_sentences = csdata.get_sentence_strings()
    print(csdata)
