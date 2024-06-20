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
        adata = sc.read_csv(HERE / 'small_data.csv').T
        csdata = cs.CSData.from_adata(adata)
        # add additional assertions to validate correct csdata creation
        assert 'CSData' in (str(csdata) + '')
