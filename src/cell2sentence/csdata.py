"""
Main data wrapper class definition
"""

#
# @author Rahul Dhodapkar
#

import numpy as np
import pandas as pd
from tqdm import tqdm

from .utils import example_function

class CSData():
    """
    Wrapper class to abstract different types of input data that can be passed
    in cell2sentence based workflows.
    """

    def __init__(self, data_path, data_path_format='arrow'):
        """
        Core constructor, CSData class contains a data path and format,
        and may also handle some buffering and data selection options.
        """
        self.data_path = data_path # path to data file / arrow format
        self.data_path_format = data_path_format # support plaintext and arrow

    def __str__(self):
        """
        Summarize CSData object as string for debugging and logging.
        """
        return "CSData Object; Path={}, Format={}".format(
            self.data_path,
            self.data_path_format)

    def to_plain_text(self, output_path):
        """
        Print data represented by CSData to a plain text file
        Arguments:
            output_path: a string representing the path to which the output file
                         should be written.
        """
        return None

    @classmethod
    def from_adata(cls, adata):
        """
        Create new CSData object from an anndata object
        """
        return cls(data_path='')


