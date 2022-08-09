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

Files used in prototyping:
- `Chick_retina_atlas_expression_matrix.csv`
- `Chick_retina_atlas_meta.csv`

## GEO Datasets
### Mouse Developmental Retina
[Manuscript Link](https://pubmed.ncbi.nlm.nih.gov/31128945/)
[Download Data Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118614)

Files used in prototyping:
- `GSE118614_10x_aggregate.mtx.gz`
- `GSE118614_barcodes.tsv.gz`
- `GSE118614_genes.tsv.gz`

### Zebrafish Developmental Retina
[Manuscript Link](https://pubmed.ncbi.nlm.nih.gov/32467236/)
[Download Data Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122680)

Files used in prototyping:
- `GSE122680_RAW.tar`


