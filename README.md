# smsk: A Snakemake skeleton to jumpstart projects

[![Build Status](https://travis-ci.org/jlanga/smsk_trinotate.svg?branch=master)](https://travis-ci.org/jlanga/smsk_trinotate)

## 1. Description

This is a workflow based on [smsk](https://github.com/jlanga/smsk) to perform the [Transdecoder](https://transdecoder.github.io/)/[Trinotate](https://trinotate.github.io/) workflow automaticaly.

## 2. First steps

Follow the contents of the `.travis.yml` file:
bash .travis/travis_before_install.sh
  - export PATH="/home/travis/miniconda3/bin:$PATH"

install:
  - bash src/install/conda_env.sh

script:
  - source activate smsk_trinotate
  - travis_wait 30 snakemake -j
1. Clone this repo

    ```sh
    git clone https://github.com/jlanga/smsk_trinotate
    cd smsk_trinotate
    ```

2. (Optional) Install miniconda and export it to your path
    
    ```sh
    bash .travis/travis_before_install.sh
    export PATH="$HOME/miniconda3/bin:$PATH"
    ```

3. Install requirements
    ```sh
    bash src/install/conda_env.sh
    ```

4. Activate the environment:

    ```sh
    source activate smsk_trinotate
    ```

4. Execute the pipeline with test data:

    ```sh
    snakemake -j
    ```

## 3. Analyzing your data

Just paste the path of your transcriptome assembly in the second line (`transcript.fa`)

Also, to increase the speed of the analysis, set the number of chunks to process to 100. To accelerate the blastp, blastx and pfam analysis, the transcriptome and proteome will be split into that number of chunks. This takes advantage of cluster analysis and also the way that blast and hmmscan paralelize.


## 3. File organization

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



## 4. Notes

- Because RNAMMER, TmHMM and SignalP require a registrations, I do not provide rules to perform those analyses.

## Bibliography

