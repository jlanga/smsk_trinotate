#!/usr/bin/env bash

# Install miniconda

if [[ -d "$HOME"/miniconda3_"$TRAVIS_OS_NAME"/bin ]]; then
    echo "miniconda already installed."
else
    echo "Installing miniconda."
    mkdir -p "$HOME/download"
    if [[ "$TRAVIS_OS_NAME" == 'osx' ]]; then
        url="https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
    else
        url="https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    fi
    wget \
        --continue \
        --output-document "$HOME"/download/miniconda_"$TRAVIS_OS_NAME".sh \
        "$url"
    chmod +x "$HOME"/download/miniconda_"$TRAVIS_OS_NAME".sh
    "$HOME"/download/miniconda_"$TRAVIS_OS_NAME".sh \
        -u \
        -b \
        -p "$HOME"/miniconda3_"$TRAVIS_OS_NAME"
    "$HOME"/miniconda3_"$TRAVIS_OS_NAME"/bin/conda clean --all --yes
fi
