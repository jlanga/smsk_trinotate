rule unzip:
    input:
        "{file}.hmm.gz"
    output:
        "{file}.hmm"
    shell:
        "gzip --decompress --keep {input}"

rule faidx_fa:
    input:
        "{file}.fa"
    output:
        "{file}.fa.fai"
    shell:
        "samtools faidx {input}"


rule faidx_fasta:
    input:
        "{file}.fasta"
    output:
        "{file}.fasta.fai"
    shell:
        "samtools faidx {input}"


rule faidx_pep:
    input:
        "{file}.pep"
    output:
        "{file}.pep.fai"
    shell:
        "samtools faidx {input}"


rule faidx_cds:
    input:
        "{file}.cds"
    output:
        "{file}.cds.fai"
    shell:
        "samtools faidx {input}"
