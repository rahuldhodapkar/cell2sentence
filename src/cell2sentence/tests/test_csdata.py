#!/usr/bin/env python
#
# Test data loading and anndata object processing
#

import anndata as ad
import scanpy as sc
from pathlib import Path
import cell2sentence as cs
import numpy as np
import random
import math

import pytest

HERE = Path(__file__).parent


class TestDataReading:
    def test_read_adata(self):
        adata = sc.read_csv(HERE / 'small_data.csv').T
        assert adata.shape == (5, 3)


class TestCoreWorkflow:
    def test_adata_to_csdata_generation(self):
        data_save_path = "/home/sr2464/palmer_scratch/C2S_Files_Syed/c2s_api_testing"

        adata = sc.read_csv(HERE / 'small_data.csv').T
        csdata = cs.CSData.from_adata(
            adata, 
            save_path=data_save_path,
            data_path_format="arrow",
            delimiter=" "
        )
        cell_sentences = csdata.get_sentence_strings()
        
        assert 'CSData' in (str(csdata) + '')
        assert csdata.feature_names == ['G1', 'G2', 'G3']
        assert csdata.data_path == data_save_path
        assert cell_sentences == {'train': ['G3', 'G1 G2', 'G1 G2 G3'], 'validation': ['G2'], 'test': ['G1 G3']}
