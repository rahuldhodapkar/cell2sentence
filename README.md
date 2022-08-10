# cell2sentence

Reframing cells as sentences of genes, ordered by expression.

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

## Getting Test Data

Test data for `cell2sentence` can be downloaded from a number of online
sequence repositories. We have 

## Broad Single Cell Portal
### Human Retina Atlas
[Manuscript Link](https://pubmed.ncbi.nlm.nih.gov/32555229/)
[Download Data Link](https://singlecell.broadinstitute.org/single_cell/study/SCP839/cell-atlas-of-the-human-fovea-and-peripheral-retina#study-download)
[Data Summary Link](https://singlecell.broadinstitute.org/single_cell/study/SCP839/cell-atlas-of-the-human-fovea-and-peripheral-retina#study-summary)
[GEO Data Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148077)

Files used in prototyping:
- `GSE148077_count_mat_donor_H1.csv.gz`
- `GSE148077_count_mat_donor_H11.csv.gz`
- `GSE148077_count_mat_donor_H2.csv.gz`
- `GSE148077_count_mat_donor_H3.csv.gz`
- `GSE148077_count_mat_donor_H4.csv.gz`
- `GSE148077_count_mat_donor_H5.csv.gz`
- `GSE148077_count_mat_donor_H9.csv.gz`

### Chick Retina Atlas
[Manuscript Link](https://pubmed.ncbi.nlm.nih.gov/33393903/)
[Download Data Link](https://singlecell.broadinstitute.org/single_cell/study/SCP1159/a-cell-atlas-of-the-chick-retina-based-on-single-cell-transcriptomics#study-download)
[Data Summary Link](https://singlecell.broadinstitute.org/single_cell/study/SCP1159/a-cell-atlas-of-the-chick-retina-based-on-single-cell-transcriptomics#study-summary)
[GEO Data Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159107)

Files used in prototyping:
- `GSE159107_E12chick_count.matrix.csv.gz`
- `GSE159107_E16chick_count.matrix.csv.gz`
- `GSE159107_E18chick_count.matrix.csv.gz`

## GEO Datasets
### Mouse Developmental Retina
[Manuscript Link](https://pubmed.ncbi.nlm.nih.gov/31128945/)
[Download Data Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118614)

Files used in prototyping:
- `GSE118614_10x_aggregate.mtx.gz`
- `GSE118614_barcodes.tsv.gz`
- `GSE118614_genes.tsv.gz`

NOTE that the genes and barcodes files available were not in the correct
formats, and needed to be cleaned with the following commands.  First move
original files,

     mv GSE118614_genes.tsv GSE118614_genes.tsv.bak
     mv GSE118614_barcodes.tsv GSE118614_barcodes.tsv.bak
     mv GSE118614_matrix.mtx GSE118614_matrix.mtx.bak

And then clean:

     tail -n +2 GSE118614_genes.tsv.bak | cut -w -f2-3 | sed s/\"//g > GSE118614_genes.tsv
     tail -n +2 GSE118614_barcodes.tsv.bak | cut -w -f2 | sed s/\"//g > GSE118614_barcodes.tsv

Additionally, the `mtx` file has been saved in the incorrect format, so
we will need to load it separately, transpose it, and save it back.

    #!/usr/bin/env python
    from scipy.io import mmread, mmwrite
    X = mmread('GSE118614_matrix.mtx.bak')
    mmwrite('GSE118614_matrix.mtx', X.T)

This mtx file is quite large, and so this will take some time.


The `.bak` files can be removed or can remain within the directory. They
will be ignored by `scanpy.read_10x_mtx`.  The `mtx`

### Zebrafish Developmental Retina
[Manuscript Link](https://pubmed.ncbi.nlm.nih.gov/32467236/)
[Download Data Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122680)

Files used in prototyping:
- `GSE122680_RAW.tar`

NOTE that after untar'ing this file, it contains several individual sequence
files, which are used separately for corpus generation 
