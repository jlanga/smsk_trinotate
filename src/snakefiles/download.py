rule download_uniprot_sprot:
    output:
        download + "uniprot_sprot.dat.gz"
    params:
        url = features["swissprot"]
    log:
        download + "uniprot_sprot.log"
    benchmark:
        download + "uniprot_sprot.json"
    conda:
        "download.yml"
    shell:
        """
        wget \
            --continue \
            --output-document {output} \
            {params.url} \
        2> {log}
        """


rule download_nog_annotations:
    output:
        download + "NOG.annotations.tsv.gz"
    params:
        url = features["NOG.annotations"]
    log:
        download + "NOG.annotations.log"
    benchmark:
        download + "NOG.annotations.json"
    conda:
        "download.yml"
    shell:
        """
        wget \
            --continue \
            --output-document {output} \
            {params.url} \
        2> {log}
        """


rule download_obo:
    output:
        download + "go-basic.obo"
    params:
        url = features["obo"]
    log:
        download + "obo.log"
    benchmark:
        download + "obo.json"
    conda:
        "download.yml"
    shell:
        """
        wget \
            --continue \
            --output-document {output} \
            {params.url} \
        2> {log}
        """


rule download_pfama:
    output:
        download + "Pfam-A.hmm.gz"
    params:
        url = features["Pfam-A"]
    log:
        download + "pfama.log"
    benchmark:
        download + "pfama.json"
    conda:
        "download.yml"
    shell:
        """
        wget \
            --continue \
            --output-document {output} \
            {params.url} \
        2> {log}
        """
