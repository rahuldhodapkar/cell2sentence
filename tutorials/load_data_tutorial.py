#!/usr/bin/env python
## load_data_tutorial.py
#
# Tutorial for loading data and generating cell sentences.
#

import numpy as np
import cell2sentence as cs
import scanpy as sc
import pandas as pd

from tqdm import tqdm
import os

os.makedirs('./calc', exist_ok=True)

# standard cell2sentence workflow:

# (1) load count matrix into anndata object depending on input format.
chick_adata = sc.read_csv('./data/chick/GSE159107_E12chick_count.matrix.csv.gz')
chick_adata = chick_adata.T # ensure that obs = cells, vars = genes.

# (2) generate sentences from anndata object.
chick_sentences = cs.transforms.generate_sentences(chick_adata)

# (3) use anndata object to build vocabulary dictionary.
chick_vocab = cs.transforms.generate_vocabulary(chick_adata)

# (4) split sentences into train, test, and validation segments.
train_chick, test_chick, val_chick = (
    cs.transforms.train_test_validation_split(chick_sentences))

