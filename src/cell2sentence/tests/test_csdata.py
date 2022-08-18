#!/usr/bin/env python
#
# Test data loading and anndata object processing
#

import anndata as ad
import scanpy as sc
from pathlib import Path
import cell2sentence as cs
import numpy as np

import pytest

HERE = Path(__file__).parent


def assert_vocab_correct(csdata):
    enc_to_feat = list(csdata.vocab.keys())
    countwords = np.repeat(0, len(enc_to_feat))
    for i in range(len(csdata.sentences)):
        for j in range(len(csdata.sentences[i])):
            countwords[csdata.sentences[i][j]] += 1

    for i, k in enumerate(csdata.vocab):
        assert countwords[i] == csdata.vocab[k]


class TestDataReading:
    def test_read_adata(self):
        adata = sc.read_csv(HERE / 'small_data.csv').T
        assert adata.shape == (5, 3)

    def test_csdata_from_adata(self):
        adata = sc.read_csv(HERE / 'small_data.csv').T
        csdata = cs.transforms.csdata_from_adata(adata)
        assert(list(csdata.vocab.keys()) == ['g1', 'g2', 'g3'])
        assert_vocab_correct(csdata)

    def test_merge_csdata(self):
        adata1 = sc.read_csv(HERE / 'small_data.csv').T
        csdata1 = cs.transforms.csdata_from_adata(adata1)

        adata2 = sc.read_csv(HERE / 'small_data_diffgenes.csv').T
        csdata2 = cs.transforms.csdata_from_adata(adata2)

        merged_csdata = cs.transforms.merge_csdata([csdata1, csdata2])
        assert_vocab_correct(merged_csdata)
