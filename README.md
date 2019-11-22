# smsk_trinotate: A Snakemake workflow to anotate a transcriptome

[![Build Status](https://travis-ci.org/jlanga/smsk_trinotate.svg?branch=master)](https://travis-ci.org/jlanga/smsk_trinotate)

## 1. Description

This is a workflow based on [smsk](https://github.com/jlanga/smsk) to perform the [Transdecoder](https://transdecoder.github.io/)/[Trinotate](https://trinotate.github.io/) workflow automaticaly.

## 2. First steps

0. Make sure you have conda and snakemake installed

1. Clone this repo

    ```sh
    git clone https://github.com/jlanga/smsk_trinotate
    cd smsk_trinotate
    ```

2. Execute the pipeline with test data:

    ```sh
    snakemake --use-conda -j 4
    ```

## 3. Analyzing your data

Just paste the path of your transcriptome assembly in the second line 
(`transcript.fa`), and the gene-transcript mapping file (`transcript.g2t.tsv`).

## 4. File organization

The hierarchy of the folder is the one described in [Good enough practices in scientific computing](https://swcarpentry.github.io/good-enough-practices-in-scientific-computing/):

The report you want will be in  `results/trinotate/trinotate.tsv`

```
smsk
├── bin: external scripts/binaries
├── data: raw data or links to backup data.
    ├── download: downloaded files required for trinotate and transdecoder
    └── db: processed databases for blast and hmmscan
├── doc: documentation.
├── README.md
├── results: processed data.
    ├── raw: links and auxiliary files
    ├── transdecoder: predicted cds and protein sequences, results from blast and hmmscan
    └── trinotate: final report, sqlite database, blastx, blastp and hmmscan results.
├── Snakefile: driver script of the project. Mostly links to src/snakefiles.
└── src: project's source code, config.yaml, snakefiles tarballs, etc.
```



## 5. Notes

- For performance reasons, `diamond` substitutes `blastp` and `blastx` both in Transdecoder and Trinotate.

- Because RNAMMER, TmHMM and SignalP require a registrations, I do not provide rules to perform those analyses.

__Update__ (2017-07-20): Rules for Rnammer, TmHMM and SignalP are provided but commented in the snakefiles. To run those analysis, uncomment rules `trinotate_signalp`, `trinotate_rnammer`, and `trinotate_tmhmm`, and uncomment also the blocks in the rule `trinotate_load` from file `src/snakefiles/trinotate.py`. Finally, modify the `config.yaml` file with the path to the `rnammer` executable.

## Bibliography
