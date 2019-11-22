# pylint:disable=syntax-error

import pandas as pd
import yaml

from snakemake.utils import min_version

min_version("5.3")

shell.prefix("set -euo pipefail;")

params = yaml.safe_load(open("params.yml", "r"))
features = yaml.safe_load(open("features.yml", "r"))
samples = pd.read_csv("samples.tsv", sep="\t")

MAX_THREADS = params["max_threads"]

singularity: "docker://continuumio/miniconda3:4.4.10"

ASSEMBLY = samples["assembly"][0]
snakefiles = "src/snakefiles/"

include: snakefiles + "generic.smk"
include: snakefiles + "folders.smk"
include: snakefiles + "clean.smk"
include: snakefiles + "raw.smk"
include: snakefiles + "download.smk"
include: snakefiles + "db.smk"
include: snakefiles + "transdecoder.smk"
include: snakefiles + "trinotate.smk"

rule all:
    input:
        # rules.raw.input,
        # rules.download.input,
        # rules.db.input,
        # rules.transdecoder.input,
        rules.trinotate.input
