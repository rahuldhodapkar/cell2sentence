"""
Functions for end-user interaction; prebuilts
"""

#
# @author Rahul Dhodapkar
#

from .csdata import CSData
from .csmodel import CSModel


def summarize_data(model, adata):
    """
    Function to generate direct insights from data. Takes a model and dataset
    and summarizes

    Arguments:
        model: a CSModel object that is fine-tuned to produce cell summaries
               from cell sentences
        adata: an AnnData object containing the data to be summarized.
    Return:
        A CSModel instance whose generate() function generates insights from
        cell sentences
    """
    return None


def generate_new_cells(model, adata, num_cells):
    """
    Generates new cells that are plausible for the provided source dataset.

    Arguments:
        model: a CSModel object that produces new cell sentences
        adata: an AnnData object containing the source dataset
    Return:
        An AnnData object containing the new cells generated.
    """
    return None


def generate_cell_perturbation_model(model, adata, perturb_df):
    """
    Generates new cells that are plausible for the provided source dataset.

    Arguments:
        model: a CSModel object to be fine-tuned.
        adata: an AnnData object containing the data to be summarized.
        perturb_df: a DataFrame object containing a mapping of perturbations
                    per cell in the adata provided.
    Return:
        A CSModel object which can be used to generate new perturbations
    """
    return None


def perturb_cells(model, adata, perturb_df):
    """
    Given a model trained to perturb new cells, generate perturbations of cells
    provided in the source adata file.

    Arguments:
        model: a CSModel object that is fine-tuned to perturb cells with a
               given set of perturbations.
        adata: an AnnData object containing the source cells to be perturbed
               in order aligned with `perturb_df`.
        perturb_df: a DataFrame object containing a mapping of perturbations
                    per cell to be applied
    Return:
        an AnnData object containing the perturbed cells.
    """
    return None
