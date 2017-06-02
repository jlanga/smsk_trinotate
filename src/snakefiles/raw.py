rule raw_link_assembly:
    input:
        config["assembly"]
    output:
        raw + "assembly.fasta"
    shell:
        "ln --symbolic "
            "$(readlink -f {input}) "
            "{output}"


rule raw_gene_to_trans_map:
    input:
        fasta = raw + "assembly.fasta"
    output:
        tsv = raw + "gene_to_trans_map.tsv"
    log:
        raw + "gene_to_trans_map.log"
    benchmark:
        raw + "gene_to_trans_map.json"
    shell:
        "get_Trinity_gene_to_trans_map.pl "
            "< {input.fasta} "
            "> {output.tsv} "
        "2> {log}"
