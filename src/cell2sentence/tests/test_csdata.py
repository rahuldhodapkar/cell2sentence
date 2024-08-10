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


# TODO: double check that we can distribute 10 cells from immune tissue dataset in GitHub repo
class TestCoreWorkflowOnImmuneTissueDataSubset:
    def setup_method(self):
        # Read in adata example object containing 10 cells from immune tissue dataset:
        # Domínguez Conde, C., et al. "Cross-tissue immune cell analysis reveals tissue-specific 
        #  features in humans." Science 376.6594 (2022): ea
        adata = sc.read_h5ad(HERE / 'immune_tissue_10cells.h5ad')
        self.save_dir = "/home/sr2464/palmer_scratch/C2S_Files_Syed/c2s_api_testing"
        self.save_name = "immune_tissue_10cells_csdata_arrow"
        
        # Define columns of adata.obs which we want to keep in cell sentence dataset
        self.adata_obs_cols_to_keep = ["cell_type", "tissue", "batch_condition", "organism"]
        
        # Create CSData object
        self.csdata = cs.CSData.from_adata(
            adata, 
            save_dir=self.save_dir,
            save_name=self.save_name,
            dataset_backend="arrow",
            sentence_delimiter=" ",
            label_col_names=self.adata_obs_cols_to_keep
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
    
    def test_arrow_dataset_splits_have_correct_number_of_samples(self):
        ds_dict = load_from_disk(self.csdata.data_path)
        assert ds_dict["train"].num_rows == 8
        assert ds_dict["validation"].num_rows == 1
        assert ds_dict["test"].num_rows == 1
    
    def test_arrow_dataset_has_correct_column_names(self):
        ds_dict = load_from_disk(self.csdata.data_path)
        train_ds = ds_dict["train"]
        assert train_ds.column_names == ["cell_name", "cell_sentence"] + self.adata_obs_cols_to_keep

    def test_arrow_dataset_saved_cell_types_correctly(self):
        ground_truth_cell_type_list_alphabetical = [
            'CD16-positive, CD56-dim natural killer cell, human',
            'CD4-positive helper T cell',
            'CD8-positive, alpha-beta memory T cell',
            'CD8-positive, alpha-beta memory T cell',
            'CD8-positive, alpha-beta memory T cell, CD45RO-positive',
            'alpha-beta T cell',
            'dendritic cell, human',
            'effector memory CD4-positive, alpha-beta T cell',
            'naive thymus-derived CD4-positive, alpha-beta T cell',
            'plasma cell'
        ]

        ds_dict = load_from_disk(self.csdata.data_path)
        cell_types = list(ds_dict["train"]["cell_type"]) + list(ds_dict["validation"]["cell_type"]) + \
            list(ds_dict["test"]["cell_type"])
        cell_types.sort()  # sorting because not sure which cells went to which splits
        assert cell_types == ground_truth_cell_type_list_alphabetical

    def test_feature_names_are_correct(self):
        # Test that ordering of feature names is correct
        assert self.csdata.feature_names[0] == "MIR1302-2HG"
        assert self.csdata.feature_names[1] == "FAM138A"
        assert self.csdata.feature_names[9] == "RP4-669L17"
    
    def test_number_of_feature_names_is_correct(self):
        assert len(self.csdata.feature_names) == 36503
    
    def test_first_train_cell_sentence(self):
        first_train_cell_sentence = self.cell_sentences["train"][0]
        first_train_cell_sentence_split = first_train_cell_sentence.split(" ")
        assert first_train_cell_sentence_split[0] == "MALAT1"
        assert first_train_cell_sentence_split[1] == "MT-ATP6"
        assert first_train_cell_sentence_split[2] == "MT-CO2"
