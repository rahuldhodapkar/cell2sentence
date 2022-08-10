#!/usr/bin/env python
#
# Serialize data for use with other external systems
#
# @author Rahul Dhodapkar <rahul.dhodapkar@yale.edu>
#

import os
import sys
from tqdm import tqdm



def xlm_prepare_outpath(outpath, species_tag, 
                        vocab,
                        train_sentences, test_sentences, val_sentences):
    """
    Write formatted data to the outpath file location, for direct processing
    by the XLM monolinguistic translation model. If creating an outpath for
    multiple species, use the same `outpath` with different `species_tag`
    values. They will not conflict so long as species_tags are appropriately
    assigned.

    Arguments:
        outpath: directory to write files to. Will create this directory
                 if it does not already exist.
        species_tag: a short string to be used as the species name in XLM.
                     Fulfills functions analaglous to language tags such as 
                     'en', 'es', or 'zh'.
        vocab: a dictionary with the vocabulary and word frequencies for 
               the species
        train_sentences: sentences to train XLM on.
        test_sentences: sentences to test XLM on.
        val_sentences: sentences to validate XLM on.

    Return:
        None
    """
    
    os.makedirs(outpath, exist_ok=True)

    print("INFO: Writing Vocabulary File", file=sys.stderr)
    vocab_fn = "{}/vocab_{}".format(outpath, species_tag)
    with open(vocab_fn, 'w') as f:
        for k in tqdm(vocab.keys()):
            print("{} {}".format(k, vocab[k]), file=f)

    print("INFO: Writing Training Sentences", file=sys.stderr)
    train_fn = "{}/train.{}".format(outpath, species_tag)
    with open(train_fn, 'w') as f:
        for l in tqdm(train_sentences):
            print(l, file=f)

    print("INFO: Writing Testing Sentences", file=sys.stderr)
    test_fn = "{}/test.{}".format(outpath, species_tag)
    with open(test_fn, 'w') as f:
        for l in tqdm(test_sentences):
            print(l, file=f)

    print("INFO: Writing Validation Sentences", file=sys.stderr)
    val_fn = "{}/valid.{}".format(outpath, species_tag)
    with open(val_fn, 'w') as f:
        for l in tqdm(val_sentences):
            print(l, file=f)
