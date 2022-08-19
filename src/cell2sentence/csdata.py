"""
Main data wrapper class definition
"""

#
# @author Rahul Dhodapkar
#

import jellyfish
import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm
import zlib


def zlib_ncd(s1, s2):
    """
    Return the zlib normalized compression distance between two strings
    """
    bs1 = bytes(s1, 'utf-8')
    bs2 = bytes(s2, 'utf-8')

    comp_cat = zlib.compress(bs1 + bs2)
    comp_bs1 = zlib.compress(bs1)
    comp_bs2 = zlib.compress(bs2)

    return ((len(comp_cat) - min(len(comp_bs1), len(comp_bs2)))
            / max(len(comp_bs1), len(comp_bs2)))

class CSData():
    """
    Lightweight wrapper class to wrap cell2sentence results.
    """

    def __init__(self, vocab, sentences, cell_names, feature_names):
        self.vocab = vocab
        self.sentences = sentences
        self.cell_names = cell_names
        self.feature_names = feature_names

    def distance_matrix(self, dist_type='zlib_ncd', prefix_len=100):
        """
        Calculate the distance matrix for the CSData object with the specified
        edit distance method. Currently supported: ("levenshtein").

        Distance caculated as d = 1 / (1 + x) where x is the similarity score.
        """

        dist_funcs = {
            'levenshtein': jellyfish.levenshtein_distance,
            'damerau_levenshtein': jellyfish.damerau_levenshtein_distance,
            'jaro': lambda x, y: 1 - jellyfish.jaro_similarity(x, y),
            'jaro_winkler':
                lambda x, y: 1 - jellyfish.jaro_winkler_similarity(x, y, longtolerance=True),
            'zlib_ncd': zlib_ncd
        }

        is_symmetric = {
            'levenshtein': True,
            'damerau_levenshtein': True,
            'jaro': True,
            'jaro_winkler': True,
            'zlib_ncd': False
        }

        mat = np.zeros(shape=(len(self.sentences), len(self.sentences)))

        for i, s_i in enumerate(tqdm(self.sentences)):
            for j, s_j in enumerate(self.sentences):
                if j < i and is_symmetric[dist_type]:
                    mat[i, j] = mat[j, i]
                    continue

                mat[i, j] = dist_funcs[dist_type](
                    s_i[:prefix_len], s_j[:prefix_len])

        return mat

    def differential_rank(self, sentence_ixs_1, sentence_ixs_2=None):
        """
        Perform differential feature rank testing given a set of sentence indexes.
        If only one group is given, the remaining sentences are automatically used
        as the comparator group.
        """

        stats_results = []
        for i, f in enumerate(self.feature_names):
            # test feature f.

            if sentence_ixs_2 is None:
                sentence_ixs_2 = set(range(len(self.sentences))).difference(
                    set(sentence_ixs_1))

            ranks_group_1 = []
            for s_ix in sentence_ixs_1:
                ranks = np.argwhere(self.sentences[s_ix] == i)
                if len(ranks) == 0:
                    ranks_group_1.append(
                        (len(self.feature_names) - len(self.sentences[s_ix])) / 2)
                else:
                    ranks_group_1.append(np.mean(ranks))

            ranks_group_2 = []
            for s_ix in sentence_ixs_2:
                ranks = np.argwhere(self.sentences[s_ix] == i)
                if len(ranks) == 0:
                    ranks_group_2.append(
                        (len(self.feature_names) - len(self.sentences[s_ix])) / 2)
                else:
                    ranks_group_2.append(np.mean(ranks))

            wilcox_stat, pval = stats.ranksums(
                x=ranks_group_1, y=ranks_group_2
            )
            stats_results.append({
                'feature': f,
                'w_stat': wilcox_stat,
                'p_val': pval,
                'mean_rank_group_1': np.mean(ranks_group_1),
                'mean_rank_group_2': np.mean(ranks_group_2)
            })
        return pd.DataFrame(stats_results)

    def generate_sentence_strings(self, delimiter=' '):
        """
        Convert internal sentence representation (arrays of ints) to traditional
        delimited character strings for integration with text-processing utilities.
        """
        if np.any([delimiter in x for x in self.feature_names]):
            raise ValueError(
                ('feature names cannot contain sentence delimiter "{}", ' +
                 'please re-format and try again').format(delimiter))

        enc_map = list(self.vocab.keys())

        joined_sentences = []
        for s in self.sentences:
            joined_sentences.append(delimiter.join(
                [enc_map[ord(x)] for x in s]
            ))

        return np.array(joined_sentences, dtype=object)
