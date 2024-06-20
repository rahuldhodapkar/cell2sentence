# cell2sentence

![cell2sentence workflow image](c2s_overview.png)

*Under construction, additional documentation is forthcoming.*

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

You will also need to install several utilities separately for running tests
and other various utilities scaffolded in the `Makefile`.

    pip install pytest
    pip install pylint
    pip install autopep8

These utilities are useful for development purposes but will not be packaged
as part of distributions.
