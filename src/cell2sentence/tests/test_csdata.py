#!/usr/bin/env python
#
# Test data loading and anndata object processing
#

import anndata as ad
import scanpy as sc
from pathlib import Path
from datasets import load_from_disk, DatasetDict
import cell2sentence as cs
import numpy as np
import random
import math
import os

import pytest

HERE = Path(__file__).parent


class TestDataReading:
    def test_read_adata(self):
        adata = sc.read_csv(HERE / 'small_data.csv').T
        assert adata.shape == (5, 3)


class TestCoreWorkflowOnDummyAdata:
    def setup_method(self):
        # Read in dummy adata object
        adata = sc.read_csv(HERE / 'small_data.csv').T
        self.save_dir = "/home/sr2464/palmer_scratch/C2S_Files_Syed/c2s_api_testing/small_data_HF_ds"
        self.save_name = "test_csdata_arrow"
        
        # Create CSData object
        self.csdata = cs.CSData.from_adata(
            adata, 
            save_dir=self.save_dir,
            save_name=self.save_name,
            dataset_backend="arrow",
            sentence_delimiter=" "
        )
        self.cell_sentences = self.csdata.get_sentence_strings()
        
    def test_csdata_string_representation(self):
        assert 'CSData' in (str(self.csdata) + '')

    def test_dataset_was_created(self):
        assert self.csdata.data_path == os.path.join(self.save_dir, self.save_name)
        assert self.csdata.dataset_backend == "arrow"
        assert os.path.exists(self.csdata.data_path)
    
    def test_arrow_dataset_created_correctly(self):
        ds_dict = load_from_disk(self.csdata.data_path)
        assert type(ds_dict) == DatasetDict

    def test_feature_names_are_correct(self):
        assert self.csdata.feature_names == ['G1', 'G2', 'G3']
    
    def test_cell_sentences_are_correct(self):
        assert self.cell_sentences == {
            'train': ['G3', 'G1 G2', 'G1 G2 G3'],
            'validation': ['G2'],
            'test': ['G1 G3']
        }
