"""
Prompt formatting class definition.
"""

#
# @authors: Syed Rizvi
#

import json
from pathlib import Path
from datasets import DatasetDict

HERE = Path(__file__).parent
SUPPORTED_TASKS = [
    "cell_type_prediction",
    "cell_type_generation",
]


def get_cell_sentence_str(ds_sample, num_genes: int = None):
    """
    Helper function for formatting cell sentences. Returns a cell sentence string containing
    a list of space-separated gene names. Caps number of genes at 'num_genes' if not None.

    Arguments:
        ds_sample: Huggingface dataset sample, assumed to follow C2S data schema.
        num_genes: if not None, integer representing number of genes to limit cell sentence to.
    """
    full_cell_sentence = ds_sample["cell_sentence"]
    if num_genes is not None:
        full_cell_sentence_gene_list = full_cell_sentence.split(" ")
        cell_sentence_str = " ".join(full_cell_sentence_gene_list[:num_genes])
        num_genes_str = str(num_genes)
    else:
        cell_sentence_str = full_cell_sentence
        num_genes_str = ""
    return cell_sentence_str, num_genes_str


class PromptFormatter():
    """
    Wrapper class to abstract different types of input data that can be passed
    in cell2sentence based workflows.
    """

    def __init__(self, task: str, top_k_genes: int):
        """
        Core constructor: PromptFormatter loads prompts for the given task and
        handles prompt formatting given a cell sentence dataset.
        """
        assert task in SUPPORTED_TASKS, "Specified finetuning task is not yet supported."
        assert top_k_genes > 0, "'top_k_genes' must be an integer > 0"
        self.task = task
        self.top_k_genes = top_k_genes

        self.prompts_dict = {}
        if task == "cell_type_prediction":
            with open(HERE / "prompts/single_cell_cell_type_prediction_prompts.json", "r") as f:
                self.prompts_dict = json.load(f)
        elif task == "cell_type_generation":
            with open(HERE / "prompts/single_cell_cell_type_conditional_generation_prompts.json", "r") as f:
                self.prompts_dict = json.load(f)
    
    def get_keys_for_task(self):
        """
        Depending on the task, this function will tell you what keys are supposed to be formatted in the
        model input and model output.
        """
        if self.task == "cell_type_prediction":
            model_input_keys = ["num_genes", "organism", "cell_type"]
            response_keys = ["cell_sentence"]
        elif self.task == "cell_type_generation":
            model_input_keys = ["num_genes", "organism", "cell_sentence"]
            response_keys = ["cell_type"]
        
        return model_input_keys, response_keys
    
    def format_one_dataset_split(self, hf_ds):
        """
        Helper function to loop through dataset samples, format prompts
        """
        model_inputs_list = []
        responses_list = []
        response_str = prompt_templates["response"][0]  # 1 response template

        # Get keys for model input and response which will need to be formatted
        model_input_keys, response_keys = self.get_keys_for_task()

        for cell_idx in range(hf_ds.num_rows):
            # TODO: introduce guardrails for number of tokens depending on LLM context length
            sample = hf_ds[cell_idx]
            
            # Get cell sentence
            single_cell_sentence_str, num_genes_str = get_cell_sentence_str(sample, num_genes=self.top_k_genes)
            sample["cell_sentence"] = single_cell_sentence_str
            sample["num_genes"] = num_genes_str

            # Select an input prompt, format keys
            model_input_str = random.choice(prompt_templates[MODEL_INPUT_KEY])
            organism_label = sample["organism"]
            cell_type_label = sample["cell_type"]

            # Format keys in model input
            for key in model_input_keys:
                model_input_str = model_input_str.replace(
                    "{" + key + "}", 
                    sample[key]
                )
            
            # Format key in response
            for key in response_keys:
                response_str = response_str.replace(
                    "{" + key + "}", 
                    sample[key]
                )

            model_inputs_list.append(model_input_str)
            responses_list.append(response_str)

        # Create formatted Huggingface dataset
        ds_split_dict = {
            "sample_type": self.task * len(hf_ds.num_rows),
            "model_input": model_inputs_list,
            "response": responses_list,
        }
        ds = Dataset.from_dict(ds_split_dict)
        return ds
    
    def format_prompts_for_dataset_dict(self, hf_ds_dict):
        """
        Function to handle formatting prompts depending on the task. Prompt formatting is done through
        a mapping function and Huggingface's Dataset .map() utility, which can batch samples to speed
        up processing.
        """
        formatted_train_ds = self.format_one_dataset_split(hf_ds=hf_ds_dict["train"])
        formatted_validation_ds = self.format_one_dataset_split(hf_ds=hf_ds_dict["validation"])
        formatted_test_ds = self.format_one_dataset_split(hf_ds=hf_ds_dict["test"])

        formatted_hf_ds_dict = DatasetDict({
            "train": formatted_train_ds,
            "validation": formatted_validation_ds,
            "test": formatted_test_ds,
        })
        return formatted_hf_ds_dict


# Debugging
if __name__ == "__main__":
    from pathlib import Path
    import scanpy as sc
    import cell2sentence as cs
    from datasets import load_from_disk

    HERE = Path(__file__).parent
    # Read in data
    adata = sc.read_h5ad(HERE / 'tests/immune_tissue_10cells.h5ad')
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

    hf_ds_dict = load_from_disk(csdata.data_path)

    task = "cell_type_prediction"
    top_k_genes = 10  # up to 10 genes
    prompt_formatter = PromptFormatter(task=task, top_k_genes=top_k_genes)

    hf_ds_dict = prompt_formatter.format_prompts_for_dataset_dict(hf_ds_dict)
    print("Done.")
