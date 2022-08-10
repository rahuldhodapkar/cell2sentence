#!/usr/bin/env python
## load_data_tutorial.py
#
# Tutorial for loading data and generating cell sentences.
#

import numpy as np
import cell2sentence as cs

chick_values = cs.transforms.read_gbc_csv(
    './data/Chick_retina_atlas_expression_matrix.csv')
chick_values['gene_names'] = np.array(['chick_' + x for x in chick_values['gene_names']])
chick_sentences = cs.transforms.trans_expression_matrix(**chick_values)
chick_vocab = cs.transforms.create_vocabulary_dict(**chick_values)


human_values = cs.transforms.read_gbc_csv(
    './data/Human_retina_combined_all_expression_matrix.csv')
human_values['gene_names'] = np.array(['human_' + x for x in human_values['gene_names']] )

human_sentences = cs.transforms.trans_expression_matrix(**human_values)
human_vocab = cs.transforms.create_vocabulary_dict(**human_values)


os.path.getsize('./data/Chick_retina_atlas_expression_matrix.csv')