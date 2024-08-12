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

For validations, we also use:
- `Human_retina_combined_all_meta.csv`

### Chick Retina Atlas
[Manuscript Link](https://pubmed.ncbi.nlm.nih.gov/33393903/)
[Download Data Link](https://singlecell.broadinstitute.org/single_cell/study/SCP1159/a-cell-atlas-of-the-chick-retina-based-on-single-cell-transcriptomics#study-download)
[Data Summary Link](https://singlecell.broadinstitute.org/single_cell/study/SCP1159/a-cell-atlas-of-the-chick-retina-based-on-single-cell-transcriptomics#study-summary)
[GEO Data Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159107)

Files used in prototyping:
- `GSE159107_E12chick_count.matrix.csv.gz`
- `GSE159107_E16chick_count.matrix.csv.gz`
- `GSE159107_E18chick_count.matrix.csv.gz`
- `Chick_retina_atlas_meta.csv`

### Macaque Foveal and Peripheral Retina
[Manuscript Link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6424338/)
[Download Data Link](https://singlecell.broadinstitute.org/single_cell/study/SCP212/molecular-specification-of-retinal-cell-types-underlying-central-and-peripheral-vision-in-primates#study-download)

Not yet used or downloaded, but could represent an interesting close relative to compare with human for
cross species translation prototyping.

Files used in prototyping:
- `Macaque_fov_AC_expression2.txt`
- `Macaque_fov_BC_expression.txt`
- `Macaque_fov_EpiImmune_expression.txt`
- `Macaque_fov_HC_expression.txt`
- `Macaque_fov_PR_expression.txt`
- `Macaque_fov_RGC_expression.txt`

## GEO Datasets
### Mouse Developmental Retina
[Manuscript Link](https://pubmed.ncbi.nlm.nih.gov/31128945/)
[Download Data Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118614)

Files used in prototyping:
- `GSE118614_10x_aggregate.mtx.gz`
- `GSE118614_barcodes.tsv.gz`
- `GSE118614_genes.tsv.gz`
- `GSE118614_barcodes.tsv.gz`

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

# Validation Data
## Human Cell Atlas Retina
[Download Link](https://singlecell.broadinstitute.org/single_cell/study/SCP737/hca-wongadultretina-retina-cumulus#study-download)
[Manuscript Link](https://www.embopress.org/doi/full/10.15252/embj.2018100811#:~:text=The%20transcriptome%20of%20human%20neural,cells%20or%20primary%20retinal%20cells.)

Files used in prototyping:
- `output.scp.matrix.mtx`
- `output.scp.features.tsv`
- `output.scp.barcodes.tsv`

We also get the cell type annotations from the human cell atlas data portal [here](https://data.humancellatlas.org/explore/projects/8185730f-4113-40d3-9cc3-929271784c2b/project-matrices)
from where we download the `WongRetinaCelltype.csv` file containing cell barcode to celltype mappings.


