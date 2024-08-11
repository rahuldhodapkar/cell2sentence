"""
Main model wrapper class definition
"""

#
# @authors: Rahul Dhodapkar, Syed Rizvi
#

# Python built-in libraries
import os
from datetime import datetime

# Third-party libraries
from datasets import load_from_disk

# Pytorch, Huggingface imports
from transformers import AutoTokenizer, AutoModelForCausalLM, TrainingArguments

# Local imports
from cell2sentence.prompt_formatter import PromptFormatter


def tokenize_all(examples, tokenizer):
    # Build list of full input prompts: model_input + response, separated by a space
    full_inputs = []
    num_samples = len(examples["model_input"])
    for sample_idx in range(num_samples):
        full_inputs.append(examples["model_input"][sample_idx] + " " + examples["response"][sample_idx])

    # Tokenize full input prompts using LLM tokenizer -> will turn prompts into input indices for LLM
    model_inputs = tokenizer(full_inputs)
    for i in range(num_samples):
        model_inputs["input_ids"][i] += [tokenizer.eos_token_id]
        model_inputs["attention_mask"][i] += [tokenizer.eos_token_id]
    model_inputs["labels"] = model_inputs["input_ids"]
    return model_inputs


class CSModel():
    """
    Wrapper class to abstract different types of input data that can be passed
    in cell2sentence based workflows.
    """

    def __init__(self, model_path_or_path, save_dir, save_name):
        """
        Core constructor, CSModel class contains a path to a model.

        Arguments:
            model_path_or_path: either a string representing a Huggingface model if 
                want to start with a default LLM, or a path to an already-trained C2S
                model on disk if want to do inference with/finetune starting from
                an already-trained C2S model.
            save_dir: directory where model should be saved to
            save_name: name to save model under (no file extension needed)
        """
        self.model_path_or_path = model_path_or_path  # path to model to load
        self.save_dir = save_dir

        # Create save path
        if not os.path.exists(save_dir):
            os.mkdir(save_dir)
        
        save_path = os.path.join(save_dir, save_name)
        self.save_path = save_path

        # Load model, save to disk
        self.tokenizer = AutoTokenizer.from_pretrained(
            model_path_or_path
        )

        model = AutoModelForCausalLM.from_pretrained(
            model_path_or_path,
            cache_dir=os.path.join(save_dir, ".cache"),
            trust_remote_code=True
        )
        model.save_checkpoint(self.save_path)  # TODO: check this

    def __str__(self):
        """
        Summarize CSData object as string for debugging and logging.
        """
        return f"CSModel Object; Path={self.model_path}"

    def fine_tune(self, csdata, task: str, top_k_genes: int = 100):
        """
        Fine tune a model using the provided CSData object data

        Arguments:
            csdata: a CSData object to be used as input for finetuning.
                    alternatively, data can be any generator of sequential
                    text that satisfies the same functional contract as
                    a CSData object.
            task: name of finetuning task (see supported tasks in prompt_formatter.py)
            top_k_genes: number of genes to use for each cell sentence
        Return:
            None: an updated CSModel is generated in-place
        """
        # Load data from csdata object
        if csdata.data_path_format == "arrow":
            hf_ds_dict = load_from_disk(csdata.data_path)
        else:
            raise NotImplementedError("Please use arrow backend implementation for training")
        
        # # Define prompt formatter, format prompts
        prompt_formatter = PromptFormatter(task=task, top_k_genes=top_k_genes)
        formatted_hf_ds_dict = prompt_formatter.format_prompts_for_dataset_dict(hf_ds_dict)

        # Tokenize data using LLM tokenizer
        # - this function applies a lambda function to tokenize each dataset split in the DatasetDict
        hf_ds_dict = hf_ds_dict.map(
            lambda batch: tokenize_all(batch, tokenizer),
            batched=True,
            load_from_cache_file=False,
            num_proc=3,
            batch_size=1000,
        )

        # Define important parameters:
        block_size = model.config.maximum_context_length  # TODO: check this this is the maximum context length of the LLM (recommended value: 4096)

        def data_collator(examples):
            # Note: this data collator assumes we are using flash attention, which does not
            #  support padding. Samples are stacked together to form blocks of input tokens.
            # concatenate samples since no padding supported
            all_input_tokens = list(
                itertools.chain(
                    *[examples[i]["input_ids"] for i in range(len(examples))]
                )
            )
            all_label_tokens = list(
                itertools.chain(*[examples[i]["labels"] for i in range(len(examples))])
            )

            assert len(all_label_tokens) == len(
                all_input_tokens
            ), "Mismatch in input ids and labels sizes"

            num_blocks = max(1, min(
                training_args.max_num_blocks,
                len(all_label_tokens) // block_size,
            ))
            assert num_blocks > 0, "Need more samples to form a block"

            all_input_tokens = torch.tensor(all_input_tokens)[: num_blocks * training_args.block_size]
            effective_block_size = min(
                training_args.block_size,
                all_input_tokens.size(all_input_tokens.dim() - 1),
            )

            if all_input_tokens.dim() == 1:
                all_input_tokens = all_input_tokens.unsqueeze(0)
            all_input_tokens = all_input_tokens.reshape(-1, effective_block_size)

            all_label_tokens = torch.tensor(all_label_tokens)[:num_blocks * training_args.block_size]
            if all_label_tokens.dim() == 1:
                all_label_tokens = all_label_tokens.unsqueeze(0)
            all_label_tokens = all_label_tokens.reshape(-1, effective_block_size)
            all_attention_mask = torch.ones_like(all_label_tokens)

            return {
                "input_ids": all_input_tokens,
                "attention_mask": all_attention_mask,
                "labels": all_label_tokens,
            }

        # Define training arguments
        datetimestamp = datetime.now().strftime('%Y-%m-%d-%H_%M_%S')
        output_dir = os.path.join(self.save_dir, datetimestamp + f"_finetune_{task}")
        train_args = TrainingArguments(
            bf16=True,
            fp16=False,
            per_device_train_batch_size=4,
            per_device_eval_batch_size=4,
            gradient_accumulation_steps=4,
            gradient_checkpointing=False,
            learning_rate=2e-4,
            load_best_model_at_end=True,
            logging_steps=50,
            logging_strategy="steps",
            lr_scheduler_type="cosine",
            num_train_epochs=10, 
            eval_steps=50,
            evaluation_strategy="steps",
            save_steps=200,
            save_strategy="steps",
            save_total_limit=3,
            warmup_ratio=0.05,
            output_dir=output_dir
        )

        # Define Trainer
        trainer = Trainer(
            model=model,
            args=train_args,
            data_collator=data_collator,
            train_dataset=train_dataset,
            eval_dataset=val_dataset,
            tokenizer=tokenizer
        )
        trainer.train()
        # TODO: change self.save_path to best model from finetuning
        print(f"Finetuning completed. Updated model saved to disk at: {self.save_path}")

        # TODO/To think about:
        # There might be many ways for people to finetune, e.g. for generation or classification, etc.
        # They will want to finetune on their own dataset of interest.
        # Need to come up with a way to specify what task they want to do. Probably good to imagine
        #  this from the end users perspective first, how would a biologist want to use this to finetune
        #  on their own data. Might be worth making a notebook that uses CSModel to train a model, and
        #  think about abstract way biologist would want to specify their task and input/output columns

        return None

    def generate(self, n=100):
        """
        Generate new data using the model.

        Arguments:
            n: the number of tokens to generate given the model supplied.
        Return:
            Text corresponding to the number `n` of tokens requested
        """
        return None

    def generate_from_prompt(self, prompt, n=100):
        """
        Generate new data using the model, starting with a given prompt.

        Arguments:
            prompt: a textual prompt.
            n: the number of tokens to generate given the model supplied.
        Return:
            Text corresponding to the number `n` of tokens requested
        """
        return None
    
    def cell_type_prediction_inference(self):
        pass


# Debugging
if __name__ == "__main__":
    from pathlib import Path

    HERE = Path(__file__).parent

    # Read in immune tissue data
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

    # Define CSModel object
    cell_type_pred_model_path = "/home/sr2464/palmer_scratch/C2S_Files_Syed/multicell_pretraining_v2_important_models/pythia-410m-multicell_v2_2024-07-28_13-55-51_checkpoint-7600_cell_type_pred"
    csmodel = CSModel(
        model_path_or_path=cell_type_pred_model_path,
        save_dir="/home/sr2464/palmer_scratch/C2S_Files_Syed/c2s_api_testing/csmodel_testing",
        save_name="cell_type_pred_pythia_410M_1"
    )

    #--- Start training from a Pythia-410M model ---#
    # # Define CSModel object
    # csmodel = CSModel(
    #     model_path_or_path="EleutherAI/pythia-410m",
    #     save_dir="/home/sr2464/palmer_scratch/C2S_Files_Syed/c2s_api_testing/csmodel_testing",
    #     save_name="training_run_1_chkpt"
    # )
    
