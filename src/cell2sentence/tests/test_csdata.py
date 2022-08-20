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


def assert_vocab_correct(csdata):
    enc_to_feat = list(csdata.vocab.keys())
    countwords = np.repeat(0, len(enc_to_feat))
    for i in range(len(csdata.sentences)):
        for j in range(len(csdata.sentences[i])):
            countwords[ord(csdata.sentences[i][j])] += 1

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

    def test_random_tiebreaks(self):
        X = np.array([[1, 1, 0],
                      [2, 2, 3]], dtype=np.int64)
        adata = ad.AnnData(
            X=X,
            obs=dict(Obs=["c1", "c2"]),
            var=dict(Feat=["g1", "g2", "g3"]),
            dtype=X.dtype
        )
        ranks_g1 = []
        for i in range(1000):
            csdata = cs.transforms.csdata_from_adata(
                adata, random_state=random.randint(0, 100))
            ranks_g1.append(csdata.sentences[0].find(chr(0)))
        assert np.mean(ranks_g1) < 0.6
        assert np.mean(ranks_g1) > 0.4


class TestSentenceSeralization:
    def test_gen_sentence_strings(self):
        adata = sc.read_csv(HERE / 'small_data.csv').T
        csdata = cs.transforms.csdata_from_adata(adata)
        strings = csdata.create_sentence_strings()
        assert len(strings) == len(csdata.sentences)


class TestSentenceProcessing:
    def test_edit_distance(self):
        adata = sc.read_csv(HERE / 'small_data.csv').T
        csdata = cs.transforms.csdata_from_adata(adata)
        mat = csdata.create_distance_matrix(dist_type='levenshtein')
        assert mat[3, 4] > mat[3, 3]
        assert mat[3, 4] == mat[4, 3]


class TestRankExtraction:
    def test_rank_extraction(self):
        adata = sc.read_csv(HERE / 'small_data.csv').T
        csdata = cs.transforms.csdata_from_adata(adata)
        rank_vec = csdata.get_rank_data_for_feature('g3')
        assert rank_vec[0] == 1
        assert rank_vec[1] == 2
        assert math.isnan(rank_vec[2])
