# pylint:disable=syntax-error

import pandas as pd
import yaml

from snakemake.utils import min_version

min_version("5.0")

shell.prefix("set -euo pipefail;")

params = yaml.load(open("params.yml", "r"))
features = yaml.load(open("features.yml", "r"))
samples = pd.read_table("samples.tsv")

singularity: "docker://continuumio/miniconda3:4.4.10"

ASSEMBLY = samples["assembly"][0]
snakefiles = "src/snakefiles/"

include: snakefiles + "generic.py"
include: snakefiles + "folders.py"
include: snakefiles + "clean.py"
include: snakefiles + "raw.py"
include: snakefiles + "download.py"
include: snakefiles + "db.py"
include: snakefiles + "transdecoder.py"
include: snakefiles + "trinotate.py"

rule all:
    input:
        # Trinotate preparations
        # # Swissprot
        # download + "uniprot_sprot.dat.gz",
        # db + prefix + ".TaxonomyIndex",
        # db + prefix + ".UniprotIndex",
        # download +  "uniprot_sprot.pep",
        # # Eggnog
        # download + "NOG.annotations.tsv.gz",
        # db + "NOG.annotations.tsv.bulk_load",
        # # Pfam-A
        # download + "Pfam-A.hmm.gz",
        # db + "Pfam-A.hmm.gz.pfam_sqlite_bulk_load",
        # # obo
        # download + "go-basic.obo",
        # db + "go-basic.obo.tab",
        # # Db
        # db + prefix + ".sqlite",
        # db + prefix + ".loaded",
        # Transdecoder
        # transdecoder + "transdecoder.cds.fai",
        # transdecoder + "transdecoder.pep.fai"
        # trinotate + "blastx.tsv",
        # trinotate + "blastp.tsv",
        # trinotate + "hmmscan.tsv",
        # trinotate + "init.txt",
        trinotate + "trinotate.tsv"
