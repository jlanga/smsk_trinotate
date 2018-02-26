#!/usr/bin/env bash

# Create conda environment
export PATH="$HOME/miniconda3_$TRAVIS_OS_NAME/bin:$PATH"
#conda env update --file environment.yml
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda install --yes -c bioconda "snakemake=4.7.0"
conda install --file .travis/packages.yml
conda clean --all --yes
