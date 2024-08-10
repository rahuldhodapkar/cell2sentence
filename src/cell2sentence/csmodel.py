"""
Main model wrapper class definition
"""

#
# @authors: Rahul Dhodapkar, Syed Rizvi
#

from datasets import load_from_disk

SUPPORTED_TASKS = [
    "cell_type_classification",
    "cell_type_generation",
]


class CSModel():
    """
    Wrapper class to abstract different types of input data that can be passed
    in cell2sentence based workflows.
    """

    def __init__(self, model_path):
        """
        Core constructor, CSModel class contains a path to a model.
        """
        self.model_path = model_path  # path to model to load

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
            task: finetuning task (currently supported; cell_type_classification 
                    and cell_type_generation)
            top_k_genes: number of genes to use for each cell sentence
        Return:
            None: an updated CSModel is generated in-place
        """
        assert task in SUPPORTED_TASKS, "Specified finetuning task is not yet supported."

        # Load data from csdata object
        if csdata.data_path_format == "arrow":
            hf_ds_dict = load_from_disk(csdata.data_path)
        else:
            raise NotImplementedError("Please use arrow backend implementation for training")
        
        # Define formatter to format prompts
        prompt_formatter = PromptFormatter(task=task, top_k_genes=top_k_genes)
        hf_ds_dict = prompt_formatter.format_prompts(hf_ds_dict)


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
            n: the number of toekns to generate given the model supplied.
        Return:
            Text corresponding to the number `n` of tokens requested
        """
        return None
