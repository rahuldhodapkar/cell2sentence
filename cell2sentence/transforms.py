#!/usr/bin/env python
#
# Transform data from common structures used in the analysis of single-cell and
# single-nucleus RNA sequencing data, to cell sentences that can be used as
# input for natural language processing tools.
#
# @author Rahul Dhodapkar <rahul.dhodapkar@yale.edu>
#

from tqdm import tqdm
import os
import numpy as np
from scipy import sparse
import sys

def read_gbc_csv(csv_filename,delimiter=',',dtype=float):
    """
    Read a gene by cell matrix csv file to a scipy sparse matrix. Assumes
    entries are numeric

    Arguments:
        csv_filename: path to the gene by cell file
        delimiter: default = ','. Field delimiter for the csv file. Can be
                   adjusted to support tsv or other similar file formats.
        dtype: default = float. Data type for csv file entries.

    Return:
        a cell by gene scipy sparse matrix
    """
    cell_names=[]
    gene_names=[]

    nonzero_gene_ixs_lst=[]
    nonzero_cell_ixs_lst=[]
    nonzero_values_lst=[]

    with open(csv_filename, 'r') as f:
        cell_names = f.readline().rstrip('\n').split(delimiter,)[1:]

        gene_ix = 0
        pbar = tqdm(total=os.path.getsize(csv_filename),
                    unit='B',
                    unit_scale=True)
        line = f.readline()
        while line:
            toks = line.rstrip('\n').split(delimiter)
            gene_names.append(toks[0])

            exp_vals = np.array([dtype(x) for x in toks[1:]], dtype=dtype)
            nonzero_values_lst.append(exp_vals[exp_vals > 0])
            nonzero_cell_ixs_lst.append(np.where(exp_vals > 0))
            nonzero_gene_ixs_lst.append(np.repeat(gene_ix, 
                                          np.sum(exp_vals > 0)))

            pbar.n = f.tell()
            pbar.refresh()
            line = f.readline()
            gene_ix += 1
        pbar.close()

    nonzero_gene_ixs=np.concatenate(nonzero_gene_ixs_lst, axis=None, dtype=int)
    nonzero_cell_ixs=np.concatenate(nonzero_cell_ixs_lst, axis=None, dtype=int)
    nonzero_values=np.concatenate(nonzero_values_lst, axis=None, dtype=dtype)

    sparse_matrix = sparse.coo_matrix(
        arg1=(nonzero_values, (nonzero_cell_ixs, nonzero_gene_ixs)),
        shape=(len(cell_names), len(gene_names))
    )

    return({
        'cell_names': np.array(cell_names),
        'gene_names': np.array(gene_names),
        'mat': sparse.csr_matrix(sparse_matrix)
    })


def create_vocabulary_dict(mat, gene_names, **kwargs):
    """
    Create a vocabulary dictionary, where each key represents a single gene
    token and the value represents the number of non-zero cells in the provided
    count matrix.

    Arguments:
        mat: a sparse CSR cell by gene matrix of raw transcript counts.
        kwargs capture other parameters to allow easy **kwargs calls.
    Return:
        a dictionary of gene vocabulary
    """
    vocabulary = {}
    gene_sums = np.sum(mat > 0, axis=0)

    for i in range(len(gene_names)):
        vocabulary[gene_names[i]] = gene_sums[0,i]

    return(vocabulary)


def trans_expression_matrix(mat, gene_names, delimiter=" ", **kwargs):
    """
    Transform expression matrix to sentences. Sentences contain gene "words"
    denoting genes with non-zero expression. Genes are ordered from highest
    expression to lowest expression.

    Arguments:
        mat: a sparse CSR cell by gene matrix of raw transcript counts.
        delimiter: default = ' '. A token delimter for the generated sentences.
        kwargs capture other parameters to allow easy **kwargs calls.
    Return:
        a list of sentences, split by delimiter.
    """
    if not isinstance(mat, sparse.csr_matrix):
        print('WARN: mat is not a scipy.sparse.csr_matrix object, attempting coercion.')
        mat = sparse.csr_matrix(mat)

    if not isinstance(gene_names, np.ndarray):
        print('WARN: gene_names is not a numpy.ndarray object, attempting coercion.')
        gene_names = np.array(gene_names)

    sentences = []
    for i in tqdm(range(mat.shape[0])):
        cols = mat.indices[mat.indptr[i]:mat.indptr[i+1]]
        vals = mat.data[mat.indptr[i]:mat.indptr[i+1]]

        words = gene_names[cols[np.argsort(-vals)]]
        sentences.append(delimiter.join(words))

    return(sentences)
