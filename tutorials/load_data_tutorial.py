#!/usr/bin/env python
## load_data_tutorial.py
#
# Tutorial for loading data and generating cell sentences.
#

import cell2sentence as cs
import os

chick_values = cs.transforms.read_gbc_csv(
    './data/Chick_retina_atlas_expression_matrix.csv')
chick_sentences = cs.transforms.trans_expression_matrix(**chick_values)
chick_vocab = cs.transforms.create_vocabulary_dict(**chick_values)


os.path.getsize('./data/Chick_retina_atlas_expression_matrix.csv')