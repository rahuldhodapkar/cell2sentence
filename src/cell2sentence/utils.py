"""
Utility functions for cell2sentence
"""

#
# @author Rahul Dhodapkar
#

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
