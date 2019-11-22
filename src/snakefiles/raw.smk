rule raw_link_assembly:
    input: ASSEMBLY
    output: RAW + "assembly.fasta"
    conda: "raw.yml"
    shell:
        "ln --symbolic $(readlink --canonicalize {input}) {output}"


rule raw_link_gene_to_trans_map:
    input: GENE_TO_TRANS_MAP
    output: RAW + "gene_to_trans_map.tsv"
    log: RAW + "gene_to_trans_map.log"
    benchmark: RAW + "gene_to_trans_map.bmk"
    conda: "raw.yml"
    shell: "ln --symbolic $(readlink --canonicalize {input}) {output}"


rule raw:
    input:
        rules.raw_link_assembly.output,
        rules.raw_link_gene_to_trans_map.output