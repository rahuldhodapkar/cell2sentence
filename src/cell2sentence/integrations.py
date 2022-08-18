"""
Serialize data for use with other external systems
"""

#
# @author Rahul Dhodapkar <rahul.dhodapkar@yale.edu>
#

import os
import sys
from tqdm import tqdm

from .transforms import train_test_validation_split


def xlm_prepare_outpath(csdata,
                        outpath,
                        species_tag,
                        params=None):
    """
    Write formatted data to the outpath file location, for direct processing
    by the XLM monolinguistic translation model. If creating an outpath for
    multiple species, use the same `outpath` with different `species_tag`
    values. They will not conflict so long as species_tags are appropriately
    assigned.

    Note that XLM requires a dictionary sorted in order of increasing
    frequency of occurence.

    Arguments:
        csdata: a CSData object from a single species to be written.
        outpath: directory to write files to. Will create this directory
                 if it does not already exist.
        species_tag: a short string to be used as the species name in XLM.
                     Fulfills functions analaglous to language tags such as
                     'en', 'es', or 'zh'.
        params: a parameter object passed to train_test_validation_split:
    Return:
        None
    """
    if params is None:
        params = {}

    train_sentences, test_sentences, val_sentences = \
        train_test_validation_split(csdata.sentences, **params)

    os.makedirs(outpath, exist_ok=True)

    print("INFO: Writing Vocabulary File", file=sys.stderr)
    vocab_fn = "{}/vocab_{}".format(outpath, species_tag)
    with open(vocab_fn, 'w') as f:
        for k in tqdm(
            sorted(
                csdata.vocab,
                key=csdata.vocab.get,
                reverse=True)):
            if csdata.vocab[k] == 0:
                continue
            print("{} {}".format(k, csdata.vocab[k]), file=f)

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
