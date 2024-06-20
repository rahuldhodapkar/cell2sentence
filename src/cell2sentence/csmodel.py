"""
Main model wrapper class definition
"""

#
# @author Rahul Dhodapkar
#

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

    def fine_tune(self, data):
        """
        Fine tune a model using the provided CSData object data

        Arguments:
            data: a CSData object to be used as input for finetuning.
        Return:
            None: an updated CSModel is generated in-place
        """
        return None

    def generate(self, n=1):
        """
        Generate new data using the existing model.

        Arguments:
            n: the number of tokens to generate given the model supplied.
        Return:
            Text corresponding to the number `n` of tokens requested
        """
        return None
