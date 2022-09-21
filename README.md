# cell2sentence

![cell2sentence workflow image](c2s_overview.png)

Reframing cells as sentences of genes, ordered by expression. Please
read the manuscript on bioRxiv for methodological details and examples.

(https://www.biorxiv.org/content/10.1101/2022.09.18.508438)

## Stable Setup

Install `cell2sentence` from PyPI with

    pip install cell2sentence

## Convert Anndata Object to Cell Sentences

After your data is loaded into a standard AnnData `adata` object, you may 
create a cell2sentence object with:

    import cell2sentence as cs

    csdata = cs.transforms.csdata_from_adata(adata)

and generate a list of cell sentences with:

    sentences = csdata.create_sentence_lists()

A tutorial script showing how to use pretrained word vectors to analyze
the `pbmc3k` dataset used by Seurat and scanpy in their guided clustering
tutorials is available at 
[`tutorials/pbmc3k_cell_sentences.py`](tutorials/pbmc3k_cell_sentences.py)

## Training Models with Cell Sentences
The `.create_sentence_lists()` and `.create_sentence_strings()` functions
can both be used to interface with a wide variety of tools. Exact
transformations required will vary from tool to too.

### gensim
As an example, some guidance
on training a `Word2Vec` model in [`gensim`](https://pypi.org/project/gensim/)
is provided here. A tutorial
from the `gensim` team is also available
[here](https://radimrehurek.com/gensim/models/word2vec.html).

For a quickstart, once you have a `csdata` object, you can run:

    import gensim

    sentences = csdata.create_sentence_lists()
    model = gensim.models.Word2Vec(sentences=sentences,
                                   vector_size=400,
                                   window=5,
                                   min_count=1,
                                   workers=4)

The model can then be queried directly, for example, to find the
top 10 genes most similar to `'CD8B'` in the embedding, you can run:

    model.wv.most_similar['CD8B']

For more details, consult the gensim documentation.

### Further Notes

As a note, the pretrained models stored in this repository are saved
instances of gensim [`KeyedVectors`](https://radimrehurek.com/gensim/models/keyedvectors.html).

If you train any models on your own data, please submit them as a pull
request or through correspondence to rahul.dhodapkar {at} yale.edu
so others can use them!  If you prototype any new uses for cell sentences,
please reach out so it can be included here.

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

## Loading Data

All data used in the bioRxiv manuscript are publicly available, and details
are outlined in the [`DATA.md`](DATA.md) file in this repository.
