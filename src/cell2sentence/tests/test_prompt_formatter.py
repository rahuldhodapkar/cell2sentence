#!/usr/bin/env python
#
# Test prompt formatting
#

import cell2sentence as cs
from cell2sentence.prompt_formatter import PromptFormatter
from pathlib import Path
from datasets import DatasetDict, load_from_disk
import scanpy as sc
import json

import pytest

HERE = Path(__file__).parent


class TestPromptFileReading:
    def test_read_adata(self):
        with open(HERE.parent / "prompts/single_cell_cell_type_conditional_generation_prompts.json", "r") as f:
            prompts = json.load(f)

        assert type(prompts) == dict


class TestCellTypePredictionPromptFormattingOnImmuneCells:
    def setup_method(self):
        # Read in dummy adata object
        adata = sc.read_h5ad(HERE / 'immune_tissue_10cells.h5ad')
        save_dir = "/home/sr2464/palmer_scratch/C2S_Files_Syed/c2s_api_testing"
        save_name = "immune_tissue_10cells_csdata_arrow"
        
        # Define columns of adata.obs which we want to keep in cell sentence dataset
        adata_obs_cols_to_keep = ["cell_type", "tissue", "batch_condition", "organism"]
        
        # Create CSData object
        csdata = cs.CSData.from_adata(
            adata, 
            save_dir=save_dir,
            save_name=save_name,
            dataset_backend="arrow",
            sentence_delimiter=" ",
            label_col_names=adata_obs_cols_to_keep
        )

        # Load cell sentence dataset
        hf_ds_dict = load_from_disk(csdata.data_path)
        
        # Format prompts for cell type prediction
        task = "cell_type_prediction"
        top_k_genes = 10  # up to 10 genes
        prompt_formatter = PromptFormatter(task=task, top_k_genes=top_k_genes)
        self.formatted_hf_ds_dict = prompt_formatter.format_prompts_for_dataset_dict(hf_ds_dict)
        
    def test_formatted_hf_ds_dict_created_correctly(self):
        assert type(self.formatted_hf_ds_dict) == DatasetDict
        assert list(self.formatted_hf_ds_dict.keys()) == ['train', 'validation', 'test']
        assert list(self.formatted_hf_ds_dict["train"].keys()) == ['sample_type', 'model_input', 'response']
        assert self.formatted_hf_ds_dict["train"].num_rows == 8

    def test_formatted_hf_ds_contains_no_braces(self):
        found_braces = False
        for ds_split in ["train", "validation", "test"]:
            ds = self.formatted_hf_ds_dict[ds_split]
            for sample_idx in range(ds.num_rows):
                sample = ds[sample_idx]
                full_input = sample["model_input"] + " " + sample["response"]
                if ("{" in full_input) or ("}" in full_input):
                    found_braces = True
        
        assert found_braces is False

    def test_formatted_hf_ds_dict_created_correctly(self):
        assert self.formatted_hf_ds_dict["train"][0]["response"] == "alpha-beta T cell."
        assert self.formatted_hf_ds_dict["train"][0]["model_input"] == (
            "The 10 gene names below are arranged by descending expression level in a Homo sapiens "
            "cell. Determine the cell type of this cell.\nCell sentence: MALAT1 MT-ATP6 MT-CO2 MT-CO1 "
            "MT-ND4 MT-CO3 MT-CYB MT-ND3 MT-ND5 MT-ND2.\nThe cell type that these genes are most "
            "commonly linked with is:"
        )


class TestCellConditionalGenerationPromptFormattingOnImmuneCells:
    def setup_method(self):
        # Read in dummy adata object
        adata = sc.read_h5ad(HERE / 'immune_tissue_10cells.h5ad')
        save_dir = "/home/sr2464/palmer_scratch/C2S_Files_Syed/c2s_api_testing"
        save_name = "immune_tissue_10cells_csdata_arrow"
        
        # Define columns of adata.obs which we want to keep in cell sentence dataset
        adata_obs_cols_to_keep = ["cell_type", "tissue", "batch_condition", "organism"]
        
        # Create CSData object
        csdata = cs.CSData.from_adata(
            adata, 
            save_dir=save_dir,
            save_name=save_name,
            dataset_backend="arrow",
            sentence_delimiter=" ",
            label_col_names=adata_obs_cols_to_keep
        )

        # Load cell sentence dataset
        hf_ds_dict = load_from_disk(csdata.data_path)
        
        # Format prompts for cell type prediction
        task = "cell_type_generation"
        top_k_genes = 10  # up to 10 genes
        prompt_formatter = PromptFormatter(task=task, top_k_genes=top_k_genes)
        self.formatted_hf_ds_dict = prompt_formatter.format_prompts_for_dataset_dict(hf_ds_dict)
        
    def test_formatted_hf_ds_dict_created_correctly(self):
        assert type(self.formatted_hf_ds_dict) == DatasetDict
        assert list(self.formatted_hf_ds_dict.keys()) == ['train', 'validation', 'test']
        assert list(self.formatted_hf_ds_dict["train"].keys()) == ['sample_type', 'model_input', 'response']
        assert self.formatted_hf_ds_dict["train"].num_rows == 8

    def test_formatted_hf_ds_contains_no_braces(self):
        found_braces = False
        for ds_split in ["train", "validation", "test"]:
            ds = self.formatted_hf_ds_dict[ds_split]
            for sample_idx in range(ds.num_rows):
                sample = ds[sample_idx]
                full_input = sample["model_input"] + " " + sample["response"]
                if ("{" in full_input) or ("}" in full_input):
                    found_braces = True
        
        assert found_braces is False

    def test_formatted_hf_ds_dict_created_correctly(self):
        assert self.formatted_hf_ds_dict["train"][0]["response"] == "MALAT1 MT-ATP6 MT-CO2 MT-CO1 MT-ND4 MT-CO3 MT-CYB MT-ND3 MT-ND5 MT-ND2."
        assert self.formatted_hf_ds_dict["train"][0]["model_input"] == (
            "Produce a list of 10 gene names in descending order of expression which represent the "
            "expressed genes of a Homo sapiens alpha-beta T cell cell.\nCell sentence:"
        )
