#!/usr/bin/env python
## load_all_species.py
#
# More representative script to load a series of experiments for different
# species and create a combined multilingual training dataset formatted for
# direct input into XLM.
#
# @author Rahul Dhodapkar <rahul.dhodapkar@yale.edu>
#

import sys
import os
from pathlib import Path
import re

import anndata
import cell2sentence as cs
import scanpy as sc

species_tags = ['human', 'mouse', 'chick', 'zebrafish']
data_dir = './data'
outpath = './calc/xlm_outpath'

def enum_csv_samples(data_dir, tag):
    files = Path("{}/{}".format(data_dir, tag)).glob('*')
    return [f for f in files if f.is_file() and not f.stem.startswith('.')]

def enum_mtx_samples(data_dir, tag):
    files = Path("{}/{}".format(data_dir, tag)).glob('*')
    return [f for f in files if f.is_dir() and not f.stem.startswith('.')]

samples_to_process = {
    'human': enum_csv_samples(data_dir, 'human'),
    'chick': enum_csv_samples(data_dir, 'chick'),
    'mouse': enum_mtx_samples(data_dir, 'mouse'),
    'zebrafish': enum_mtx_samples(data_dir, 'zebrafish')
}

def read_csv_sample(s):
    adata = sc.read_csv(s)
    return(adata.T)

def read_mtx_sample(s):
    """
    Many 10x mtx-formatted files have an identifier prefix, such as 
    "samplename_matrix.mtx.gz". This snippet automatically detects such
    prefixes, and does nothing if no prefix is present.

    ***NOTE*** that due to a bug in scanpy, legacy 10x matrix files must
    be unzipped prior to reading, while v3 10x matrix files must be zipped
    and have the *.gz extension.
    """
    mtx_file =  list(Path(s).glob('*.mtx'))
    if (len(mtx_file) != 1):
        print("ERROR: improperly formatted 10x directory ({})".format(s))
        sys.exit(1)
    result = re.match(r'(\w*)matrix.mtx', mtx_file[0].name)
    prefix = result.group(1)
    adata = sc.read_10x_mtx(s, prefix=prefix)
    return(adata)

read_funcs = {
    'human': read_csv_sample,
    'chick': read_csv_sample,
    'mouse': read_mtx_sample,
    'zebrafish': read_mtx_sample
}

for tag in species_tags:
    print("INFO: loading data for species [{}]".format(tag))
    adata_objs = []
    for s in samples_to_process[tag]:
        adata = read_funcs[tag](s)
        adata_objs.append(adata)
    print("INFO: joining data for species [{}]".format(tag))
    adata_combined = anndata.concat(adata_objs, axis=0)
    print("INFO: generating sentences for species [{}]".format(tag))
    sentences = cs.transforms.generate_sentences(adata_combined)
    print("INFO: building vocabulary for species [{}]".format(tag))
    vocab = cs.transforms.generate_vocabulary(adata_combined)
    print("INFO: formatting for XLM integration for species [{}]".format(tag))
    train, test, val = cs.transforms.train_test_validation_split(
        sentences, train_pct=0.8, test_pct=0.1, val_pct=0.1)
    cs.integrations.xlm_prepare_outpath(
        outpath=outpath, species_tag=tag,
        vocab=vocab,
        train_sentences=train,
        test_sentences=test,
        val_sentences=val)
    print("INFO: species [{}] complete".format(tag))

print("All done!")
