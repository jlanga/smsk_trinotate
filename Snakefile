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
        # DOWNLOAD + "uniprot_sprot.dat.gz",
        # DB + prefix + ".TaxonomyIndex",
        # DB + prefix + ".UniprotIndex",
        # DOWNLOAD +  "uniprot_sprot.pep",
        # # Eggnog
        # DOWNLOAD + "NOG.annotations.tsv.gz",
        # DB + "NOG.annotations.tsv.bulk_load",
        # # Pfam-A
        # DOWNLOAD + "Pfam-A.hmm.gz",
        # DB + "Pfam-A.hmm.gz.pfam_sqlite_bulk_load",
        # # obo
        # DOWNLOAD + "go-basic.obo",
        # DB + "go-basic.obo.tab",
        # # Db
        # DB + prefix + ".sqlite",
        # DB + prefix + ".loaded",
        # Transdecoder
        # TRANSDECODER + "transdecoder.cds.fai",
        # TRANSDECODER + "transdecoder.pep.fai"
        # TRINOTATE + "blastx.tsv",
        # TRINOTATE + "blastp.tsv",
        # TRINOTATE + "hmmscan.tsv",
        # TRINOTATE + "init.txt",
        TRINOTATE + "trinotate.tsv"
