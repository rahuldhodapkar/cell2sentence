#!/usr/bin/env python
#
# Test prompt formatting
#

import cell2sentence as cs
from cell2sentence.prompt_formatter import PromptFormatter
from pathlib import Path
import json

import pytest

HERE = Path(__file__).parent


class TestPromptFileReading:
    def test_read_adata(self):
        with open(HERE / "prompts/single_cell_cell_type_conditional_generation_prompts.json", "r") as f:
            prompts = json.load(f)

        assert type(prompts) == dict


class TestPromptFormattingDummySentences:
    def setup_method(self):
        # Read in dummy adata object
        adata = sc.read_csv(HERE / 'small_data.csv').T
        self.save_dir = "/home/sr2464/palmer_scratch/C2S_Files_Syed/c2s_api_testing/small_data_HF_ds"
        self.save_name = "test_csdata_arrow"
        
        # Create CSData object
        self.csdata = cs.CSData.from_adata(
            adata, 
            save_dir=self.save_dir,
            save_name=self.save_name,
            dataset_backend="arrow",
            sentence_delimiter=" "
        )
        
    def test_cell_type_pred_prompt_formatting(self):
        hf_ds_dict = load_from_disk(self.csdata.data_path)

        task = "cell_type_prediction"
        top_k_genes = 10  # up to 10 genes
        prompt_formatter = PromptFormatter(task=task, top_k_genes=top_k_genes)
        hf_ds_dict = prompt_formatter.format_prompts_for_dataset_dict(hf_ds_dict)