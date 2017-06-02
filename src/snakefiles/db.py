rule db_parse_uniprot_sprot:
    input:
        download + "uniprot_sprot.dat.gz"
    output:
        taxonomy = db + "trinotate.TaxonomyIndex",
        uniprot = db +  "trinotate.UniprotIndex",
        pep = download + "uniprot_sprot.pep"
    params:
        pep_tmp =  download + "uniprot_sprot.dat.gz.pep",
        prefix = db + "trinotate"
    log:
        db + "parse_uniprot_sprot.log"
    benchmark:
        db + "parse_uniprot_sprot.josn"
    shell:
        "EMBL_swissprot_parser.pl "
            "{input} "
            "{params.prefix} "
        "2> {log}; "
        "mv {params.pep_tmp} {output.pep}"



rule db_obo_to_tab:
    input:
        download + "go-basic.obo"
    output:
        db + "go-basic.obo.tab"
    log:
        db + "obo_to_tab.log"
    benchmark:
        db + "obo_to_tab.json"
    shell:
        "obo_to_tab.pl {input} > {output} 2> {log}"


rule db_parse_nog:
    input:
        tab = download + "NOG.annotations.tsv.gz"
    output:
        db + "NOG.annotations.tsv.bulk_load"
    log:
        db + "parse_nog.log"
    benchmark:
        db + "parse_nog.json"
    shell:
        "(gzip -dc {input} "
        "| print.pl 1 5 "
        "> {output} ) 2> {log}"


rule db_parse_pfam:
    input:
        download + "Pfam-A.hmm.gz"
    output:
        db + "Pfam-A.hmm.gz.pfam_sqlite_bulk_load"
    params:
        tmp = download + "Pfam-A.hmm.gz.pfam_sqlite_bulk_load"
    log:
        db + "parse_pfam.log"
    benchmark:
        db + "parse_pfam.json"
    shell:
        "PFAM_dat_parser.pl {input} 2> {log}; "
        "mv {params.tmp} {output} 2>> {log}"



rule db_hmmpress_pfama:
    """
    Format the HMM Pfam-A database
    """
    input:
        hmm_gz = download + "Pfam-A.hmm.gz"
    output:
        hmm = db + "Pfam-A.hmm",
        other = expand(
            db + "Pfam-A.hmm.{extension}",
            extension = "h3i h3f h3p".split()
        )
    log:
        db + "hmmpress_pfama.log"
    benchmark:
        db + "hmmpress_pfama.json"
    shell:
        "gzip "
            "--decompress --keep --stdout "
            "{input.hmm_gz} "
        "> {output.hmm} "
        "2> {log}; "
        "hmmpress {output.hmm} 2>> {log} 1>&2"


rule db_makeblastdb_uniprot_sprot:
    """
    Make the SwissProt filtered database
    """
    input:
        download + "uniprot_sprot.pep"
    output:
        touch(db + "uniprot_sprot")
    threads:
        1
    log:
        db + "makeblastdb_uniprot_sprot.log"
    benchmark:
        db + "makeblastdb_uniprot_sprot.json"
    shell:
        "makeblastdb "
            "-dbtype prot "
            "-title {output} "
            "-out {output} "
            "-parse_seqids "
            "-in {input} "
        "2> {log} 1>&2"
