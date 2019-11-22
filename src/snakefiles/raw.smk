rule raw_link_assembly:
    input:
        ASSEMBLY
    output:
        RAW + "assembly.fasta"
    conda:
        "raw.yml"
    shell:
        """
        ln --symbolic \
            $(readlink --canonicalize {input}) \
            {output}
        """


rule raw_gene_to_trans_map:
    input:
        fasta = RAW + "assembly.fasta"
    output:
        tsv = RAW + "gene_to_trans_map.tsv"
    log:
        RAW + "gene_to_trans_map.log"
    benchmark:
        RAW + "gene_to_trans_map.bmk"
    conda:
        "raw.yml"
    shell:
        """
        get_Trinity_gene_to_trans_map.pl \
            < {input.fasta} \
            > {output.tsv} \
        2> {log}
        """


rule raw:
    input:
        rules.raw_link_assembly.output,
        rules.raw_gene_to_trans_map.output