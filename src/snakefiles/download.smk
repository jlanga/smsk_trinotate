rule download_uniprot_sprot:
    output:
        DOWNLOAD + "uniprot_sprot.dat.gz"
    params:
        url = features["swissprot"]
    log:
        DOWNLOAD + "uniprot_sprot.log"
    benchmark:
        DOWNLOAD + "uniprot_sprot.bmk"
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
        DOWNLOAD + "NOG.annotations.tsv.gz"
    params:
        url = features["NOG.annotations"]
    log:
        DOWNLOAD + "NOG.annotations.log"
    benchmark:
        DOWNLOAD + "NOG.annotations.bmk"
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
        DOWNLOAD + "go-basic.obo"
    params:
        url = features["obo"]
    log:
        DOWNLOAD + "obo.log"
    benchmark:
        DOWNLOAD + "obo.bmk"
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
        DOWNLOAD + "Pfam-A.hmm.gz"
    params:
        url = features["Pfam-A"]
    log:
        DOWNLOAD + "pfama.log"
    benchmark:
        DOWNLOAD + "pfama.bmk"
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
