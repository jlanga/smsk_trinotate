CHUNKS = params["chunks"]["transdecoder"]


rule transdecoder_longorfs:
    """
    Predict proteins according to the presence of long ORFs
    """
    input:
        fasta = RAW + "assembly.fasta",
        tsv = RAW + "gene_to_trans_map.tsv"
    output:
        "assembly.fasta.transdecoder_dir/longest_orfs.pep"
    log:
        TRANSDECODER + "longorfs.log"
    benchmark:
        TRANSDECODER + "longorfs.bmk"
    conda:
        "transdecoder.yml"
    shell:
        """
        TransDecoder.LongOrfs \
            -t {input.fasta} \
            --gene_trans_map {input.tsv} \
        2> {log} 1>&2
        """


rule transdecoder_split_longest_orfs:
    """
    Split the headers from transdecoder_longest_orfs into multiple files
    """
    input:
        fai = "assembly.fasta.transdecoder_dir/longest_orfs.pep.fai"
    output:
        expand(
            TRANSDECODER + "chunks/longest_orfs_{chunk_id}.tsv",
            chunk_id=['{0:05d}'.format(x) for x in range(0, CHUNKS)]
        )
    params:
        number_of_chunks = CHUNKS
    log:
        TRANSDECODER + "split_longest_orfs.log"
    benchmark:
        TRANSDECODER + "split_longest_orfs.bmk"
    conda:
        "transdecoder.yml"
    shell:
        """
        split \
            --number l/{params.number_of_chunks} \
            --numeric-suffixes \
            --suffix-length 5 \
            --additional-suffix .tsv \
            {input.fai} \
            {TRANSDECODER}/chunks/longest_orfs_ \
        2> {log}
        """


rule transdecoder_hmmscan_chunk:
    """
    hmmscan over one chunk
    """
    input:
        pep = "assembly.fasta.transdecoder_dir/longest_orfs.pep",
        fai = "assembly.fasta.transdecoder_dir/longest_orfs.pep.fai",
        chunk = TRANSDECODER + "chunks/longest_orfs_{chunk_id}.tsv",
        hmm = DB + "Pfam-A.hmm"
    output:
        tsv = TRANSDECODER + "hmmscan/longest_orfs_{chunk_id}.tsv"
    log:
        TRANSDECODER + "hmmscan/longest_orfs_{chunk_id}.log"
    benchmark:
        TRANSDECODER + "hmmscan/longest_orfs_{chunk_id}.bmk"
    conda:
        "transdecoder.yml"
    shell:
        """
        cut -f 1 {input.chunk} \
        | xargs samtools faidx {input.pep} \
        | hmmscan \
            --cpu {threads} \
            --domtblout {output.tsv} \
            {input.hmm} \
            - \
        2> {log} 1>&2
        """


rule transdecoder_hmmscan_merge:
    """
    Merge hmmscan results into one file
    """
    input:
        expand(
            TRANSDECODER + "hmmscan/longest_orfs_{chunk_id}.tsv",
            chunk_id=['{0:05d}'.format(x) for x in range(0, CHUNKS)]
        )
    output:
        tsv = TRANSDECODER + "hmmscan.tsv"
    log:
        TRANSDECODER + "hmmscan_merge.log"
    benchmark:
        TRANSDECODER + "hmmscan_merge.bmk"
    conda:
        "transdecoder.yml"
    shell:
        """cat {input} > {output} 2> {log}"""


rule transdecoder_blastp_chunk:
    """
    Run blastp of each chunk
    """
    input:
        pep = "assembly.fasta.transdecoder_dir/longest_orfs.pep",
        fai = "assembly.fasta.transdecoder_dir/longest_orfs.pep.fai",
        chunk = TRANSDECODER + "chunks/longest_orfs_{chunk_id}.tsv",
        blast_db = DB + "uniprot_sprot"
    output:
        tsv = TRANSDECODER + "blastp/longest_orfs_{chunk_id}.tsv"
    log:
        TRANSDECODER + "blastp/longest_orfs_{chunk_id}.log"
    benchmark:
        TRANSDECODER + "blastp/longest_orfs_{chunk_id}.bmk"
    conda:
        "transdecoder.yml"
    shell:
        """
        cut -f 1 {input.chunk} \
        | xargs samtools faidx {input.pep} \
        | blastp \
            -db {input.blast_db} \
            -max_target_seqs 1 \
            -outfmt 6 \
            -evalue 1e-5 \
            -out {output.tsv} \
        2> {log} 1>&2
        """


rule transdecoder_blastp_merge:
    """
    Merge results from the different blastps
    """
    input:
        expand(
            TRANSDECODER + "blastp/longest_orfs_{chunk_id}.tsv",
            chunk_id=[
                '{0:05d}'.format(x)
                for x in range(0, CHUNKS)
            ]
        )
    output:
        tsv = TRANSDECODER + "blastp.tsv"
    log:
        TRANSDECODER + "blastp_merge.log"
    benchmark:
        TRANSDECODER + "blastp_merge.bmk"
    conda:
        "transdecoder.yml"
    shell:
        "cat {input} > {output} 2> {log}"


rule transdecoder_predict:
    """
    Join results from blast and hmmr to predict coding sequences
    """
    input:
        fasta = RAW + "assembly.fasta",
        pfam_tsv = TRANSDECODER + "hmmscan.tsv",
        blastp_tsv = TRANSDECODER + "blastp.tsv"
    output:
        bed = TRANSDECODER + "transdecoder.bed",
        cds = TRANSDECODER + "transdecoder.cds",
        gff3 = TRANSDECODER + "transdecoder.gff3",
        pep = TRANSDECODER + "transdecoder.pep",
    params:
        dir = "assembly.fasta.transdecoder_dir",
        bed = "assembly.fasta.transdecoder.bed",
        cds = "assembly.fasta.transdecoder.cds",
        gff3 = "assembly.fasta.transdecoder.gff3",
        pep = "assembly.fasta.transdecoder.pep",
        checkpoints = "assembly.fasta.transdecoder_dir.__checkpoints"
    threads:
        24
    log:
        TRANSDECODER + "predict.log"
    benchmark:
        TRANSDECODER + "predict.bmk"
    conda:
        "transdecoder.yml"
    shell:
        """
        TransDecoder.Predict \
            -t {input.fasta} \
            --retain_pfam_hits {input.pfam_tsv} \
            --retain_blastp_hits {input.blastp_tsv} \
            --cpu {threads} \
        2> {log} 1>&2

        mv {params.bed} {output.bed} 2>> {log} 1>&2
        mv {params.cds} {output.cds} 2>> {log} 1>&2
        mv {params.gff3} {output.gff3} 2>> {log} 1>&2
        mv {params.pep} {output.pep} 2>> {log} 1>&2
        rm -rf {params.dir} {params.checkpoints} 2>> {log} 1>&2
        """
