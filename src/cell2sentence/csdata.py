#!/usr/bin/env python
#
# @author Rahul Dhodapkar
#


class CSData():
    """
    Lightweight wrapper class to wrap cell2sentence results.
    """

    def __init__(self, vocab, sentences, cell_names, feature_names):
        self.vocab = vocab
        self.sentences = sentences
        self.cell_names = cell_names
        self.feature_names = feature_names
