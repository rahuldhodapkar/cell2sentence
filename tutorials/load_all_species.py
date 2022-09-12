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
import pandas as pd
from tqdm import tqdm

import anndata
import cell2sentence as cs
import scanpy as sc

import numpy as np


#species_tags = ['human', 'mouse', 'chick', 'zebrafish', 'macaque']
species_tags = ['macaque']

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
    'zebrafish': enum_mtx_samples(data_dir, 'zebrafish'),
    'macaque': enum_csv_samples(data_dir, 'macaque')
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


def clean_adata(adata, tag):
    """
    Ensure that anndata object `adata` is correctly formatted for use by
    cell2sentence.
    """
    adata.var_names = pd.Index([tag + '_' + x.replace(' ', '_') for x in adata.var_names],dtype=object)
    return(adata)


read_funcs = {
    'human': read_csv_sample,
    'chick': read_csv_sample,
    'mouse': read_mtx_sample,
    'zebrafish': read_mtx_sample,
    'macaque': read_csv_sample
}

csdata_lst = [] # collect names to single vocabulary file
for tag in species_tags:
    print("INFO: loading data for species [{}]".format(tag))
    adata_objs = []
    for s in samples_to_process[tag]:
        adata = clean_adata(adata=read_funcs[tag](s), tag=tag)
        adata_objs.append(adata)
    print("INFO: joining data for species [{}]".format(tag))
    adata_combined = anndata.concat(adata_objs, axis=0)
    sc.pp.highly_variable_genes(adata_combined, n_top_genes=2000, flavor='seurat_v3')
    adata_hvg = adata_combined[:,adata_combined.var.highly_variable]
    print("INFO: Saving combined matrix for species [{}]".format(tag))
    save_df = adata_combined[:,adata_combined.var.highly_variable].to_df()
    save_df.to_csv('{}/{}_hvg_expression.csv'.format(outpath, tag))
    print("INFO: generating sentences for species [{}]".format(tag))
    csdata_combined = cs.transforms.csdata_from_adata(adata_hvg, prefix_len=30)
    cs.integrations.xlm_prepare_outpath(
        csdata_combined,
        outpath=outpath,
        species_tag=tag)
    csdata_lst.append(csdata_combined)
    print("INFO: species [{}] complete".format(tag))

print("INFO: saving merged vocabulary file")
csdata_merged = cs.transforms.merge_csdata(csdata_lst)
merged_vocab = csdata_merged.vocab
with open('{}/vocab'.format(outpath), 'w') as f:
    for k in tqdm(sorted(merged_vocab, key=merged_vocab.get, reverse=True)):
        if merged_vocab[k] == 0:
            continue
        print("{} {}".format(k, merged_vocab[k]), file=f)

print("All done!")
