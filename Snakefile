shell.prefix("set -euo pipefail;")
configfile: "src/config.yaml"

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
