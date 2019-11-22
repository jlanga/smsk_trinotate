rule db_parse_uniprot_sprot:
    input:
        DOWNLOAD + "uniprot_sprot.dat.gz"
    output:
        taxonomy = DB + "trinotate.TaxonomyIndex",
        uniprot = DB + "trinotate.UniprotIndex",
        pep = DOWNLOAD + "uniprot_sprot.pep"
    params:
        pep_tmp = DOWNLOAD + "uniprot_sprot.dat.gz.pep",
        prefix = DB + "trinotate"
    log:
        DB + "parse_uniprot_sprot.log"
    benchmark:
        DB + "parse_uniprot_sprot.bmk"
    conda:
        "db.yml"
    shell:
        "EMBL_swissprot_parser.pl {input} {params.prefix} "
        "2> {log}; "
        "mv {params.pep_tmp} {output.pep}"


rule db_obo_to_tab:
    input:
        DOWNLOAD + "go-basic.obo"
    output:
        DB + "go-basic.obo.tab"
    log:
        DB + "obo_to_tab.log"
    benchmark:
        DB + "obo_to_tab.bmk"
    conda:
        "db.yml"
    shell:
        "obo_to_tab.pl {input} > {output} 2> {log}"


rule db_parse_nog:
    input:
        tab = DOWNLOAD + "NOG.annotations.tsv.gz"
    output:
        DB + "NOG.annotations.tsv.bulk_load"
    log:
        DB + "parse_nog.log"
    benchmark:
        DB + "parse_nog.bmk"
    conda:
        "db.yml"
    shell:
        "(gzip -dc {input} "
        "| print.pl 1 5 "
        "> {output} ) 2> {log}"


rule db_parse_pfam:
    input:
        DOWNLOAD + "Pfam-A.hmm.gz"
    output:
        DB + "Pfam-A.hmm.gz.pfam_sqlite_bulk_load"
    params:
        tmp = DOWNLOAD + "Pfam-A.hmm.gz.pfam_sqlite_bulk_load"
    log:
        DB + "parse_pfam.log"
    benchmark:
        DB + "parse_pfam.bmk"
    conda:
        "db.yml"
    shell:
        "PFAM_dat_parser.pl {input} 2> {log}; "
        "mv {params.tmp} {output} 2>> {log}"


rule db_hmmpress_pfama:
    """
    Format the HMM Pfam-A database
    """
    input:
        hmm_gz = DOWNLOAD + "Pfam-A.hmm.gz"
    output:
        hmm = DB + "Pfam-A.hmm",
        other = expand(
            DB + "Pfam-A.hmm.{extension}",
            extension="h3i h3f h3p".split()
        )
    log:
        DB + "hmmpress_pfama.log"
    benchmark:
        DB + "hmmpress_pfama.bmk"
    conda:
        "db.yml"
    shell:
        "gzip --decompress --keep --stdout {input.hmm_gz} "
        "> {output.hmm} 2> {log}; "
        "hmmpress {output.hmm} 2>> {log} 1>&2; "
        "cat /dev/null > {output.hmm} 2>> {log} 1>&2"


rule db_diamond_makedb_uniprot_sprot:
    """
    Make the SwissProt database for Diamond
    """
    input: DOWNLOAD + "uniprot_sprot.pep"
    output: DB + "uniprot_sprot.dmnd"
    params: DB + "uniprot_sprot"
    threads: 1
    log: DB + "diamond_makedb_unisprot_sprot.log"
    benchmark: DB + "diamond_makedb_unisprot_sprot.bmk"
    conda: "db.yml"
    shell: "diamond makedb --in {input} --db {params} 2> {log} 1>&2"