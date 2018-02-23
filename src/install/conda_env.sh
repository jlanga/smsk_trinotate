#!/usr/bin/env bash

# Create conda environment
export PATH="$HOME/miniconda3_$TRAVIS_OS_NAME/bin:$PATH"
#conda env update --file environment.yml
conda install --yes -c bioconda snakemake
conda clean --all --yes
