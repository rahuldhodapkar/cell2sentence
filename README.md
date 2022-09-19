# cell2sentence

![cell2sentence workflow image](c2s_overview.png)

Reframing cells as sentences of genes, ordered by expression.

## Stable Setup

Install `cell2sentence` from PyPI with

    pip install cell2sentence


## Development Setup

Create a conda environment using `python3` using 
[anaconda](https://docs.anaconda.com/anaconda/install/) with:

    conda create -n cell2sentence python=3.8

and activate the environment with

    conda activate cell2sentence

finally, you can install the latest development version of `cell2sentence` by
running

    make install

which simply uses `pip -e`.
