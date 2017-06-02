rule transdecoder_longorfs:
    """
    Predict proteins according to the presence of long ORFs
    """
    input:
        fasta = raw + "assembly.fasta",
        tsv = raw + "gene_to_trans_map.tsv"
    output:
        "assembly.fasta.transdecoder_dir/longest_orfs.pep",
        temp(
            "assembly.fasta.transdecoder_dir/"
        )
    log:
        transdecoder + "longorfs.log"
    benchmark:
        transdecoder + "longorfs.json"
    shell:
        "TransDecoder.LongOrfs "
            "-t {input.fasta} "
            "--gene_trans_map {input.tsv} "
        "2> {log} 1>&2"



rule transdecoder_split_longest_orfs:
    """
    Split the headers from transdecoder_longest_orfs into multiple files
    """
    input:
        fai = "assembly.fasta.transdecoder_dir/longest_orfs.pep.fai"
    output:
        expand(
            transdecoder + "chunks/longest_orfs_{chunk_id}.tsv",
            chunk_id = ['{0:05d}'.format(x) for x in range(0, config["number_of_chunks"]["transdecoder"])]
        )
    params:
        number_of_chunks = config["number_of_chunks"]["transdecoder"]
    log:
        transdecoder + "split_longest_orfs.log"
    benchmark:
        transdecoder + "split_longest_orfs.json"
    shell:
        "split "
            "--number l/{params.number_of_chunks} "
            "--numeric-suffixes "
            "--suffix-length 5 "
            "--additional-suffix .tsv "
            "{input.fai} "
            "{transdecoder}/chunks/longest_orfs_ "
        "2> {log}"



rule transdecoder_hmmscan_chunk:
    """
    hmmscan over one chunk
    """
    input:
        pep = "assembly.fasta.transdecoder_dir/longest_orfs.pep",
        fai = "assembly.fasta.transdecoder_dir/longest_orfs.pep.fai",
        chunk = transdecoder + "chunks/longest_orfs_{chunk_id}.tsv",
        hmm = db + "Pfam-A.hmm"
    output:
        tsv = transdecoder + "hmmscan/longest_orfs_{chunk_id}.tsv"
    log:
        transdecoder + "hmmscan/longest_orfs_{chunk_id}.log"
    benchmark:
        transdecoder + "hmmscan/longest_orfs_{chunk_id}.json"
    shell:
        "cut -f 1 {input.chunk} "
        "| xargs samtools faidx {input.pep} "
        "| hmmscan "
            "--domtblout {output.tsv} "
            "{input.hmm} "
            "- "
        "2> {log} 1>&2"



rule transdecoder_hmmscan_merge:
    """
    Merge hmmscan results into one file
    """
    input:
        expand(
            transdecoder + "hmmscan/longest_orfs_{chunk_id}.tsv",
            chunk_id = ['{0:05d}'.format(x) for x in range(0, config["number_of_chunks"]["transdecoder"])]
        )
    output:
        tsv = transdecoder + "hmmscan.tsv"
    log:
        transdecoder + "hmmscan_merge.log"
    benchmark:
        transdecoder + "hmmscan_merge.json"
    shell:
        "cat {input} > {output} 2> {log}"



rule transdecoder_blastp_chunk:
    """
    Run blastp of each chunk
    """
    input:
        pep = "assembly.fasta.transdecoder_dir/longest_orfs.pep",
        fai = "assembly.fasta.transdecoder_dir/longest_orfs.pep.fai",
        chunk = transdecoder + "chunks/longest_orfs_{chunk_id}.tsv",
        db = db + "uniprot_sprot"
    output:
        tsv = transdecoder + "blastp/longest_orfs_{chunk_id}.tsv"
    log:
        transdecoder + "blastp/longest_orfs_{chunk_id}.log"
    benchmark:
        transdecoder + "blastp/longest_orfs_{chunk_id}.json"
    shell:
        "cut -f 1 {input.chunk} "
        "| xargs samtools faidx {input.pep} "
        "| blastp "
            "-db {input.db} "
            "-max_target_seqs 1 "
            "-outfmt 6 "
            "-evalue 1e-5 "
            "-out {output.tsv} "
        "2> {log} 1>&2"



rule transdecoder_blastp_merge:
    """
    Merge results from the different blastps
    """
    input:
        expand(
            transdecoder + "blastp/longest_orfs_{chunk_id}.tsv",
            chunk_id = ['{0:05d}'.format(x) for x in range(0, config["number_of_chunks"]["transdecoder"])]
        )
    output:
        tsv = transdecoder + "blastp.tsv"
    log:
        transdecoder + "blastp_merge.log"
    benchmark:
        transdecoder + "blastp_merge.json"
    shell:
        "cat {input} > {output} 2> {log}"



rule transdecoder_predict:
    """
    Join results from blast and hmmr to predict coding sequences
    """
    input:
        fasta = raw + "assembly.fasta",
        pfam_tsv = transdecoder + "hmmscan.tsv",
        blastp_tsv = transdecoder + "blastp.tsv",
        folder = "assembly.fasta.transdecoder_dir/"
    output:
        bed  = transdecoder + "transdecoder.bed",
        cds  = transdecoder + "transdecoder.cds",
        gff3 = transdecoder + "transdecoder.gff3",
        pep  = transdecoder + "transdecoder.pep",
    params:
        dir  = "assembly.fasta.transdecoder_dir",
        bed  = "assembly.fasta.transdecoder.bed",
        cds  = "assembly.fasta.transdecoder.cds",
        gff3 = "assembly.fasta.transdecoder.gff3",
        pep  = "assembly.fasta.transdecoder.pep",
    threads:
        24
    log:
        transdecoder + "predict.log"
    benchmark:
        transdecoder + "predict.json"
    shell:
        "TransDecoder.Predict "
            "-t {input.fasta} "
            "--retain_pfam_hits {input.pfam_tsv} "
            "--retain_blastp_hits {input.blastp_tsv} "
            "--cpu {threads} "
        "2> {log} 1>&2 ;"
        "mv {params.bed} {output.bed} 2>> {log} 1>&2; "
        "mv {params.cds} {output.cds} 2>> {log} 1>&2; "
        "mv {params.gff3} {output.gff3} 2>> {log} 1>&2; "
        "mv {params.pep} {output.pep} 2>> {log} 1>&2; "
        "rm -rf {params.dir} 2>> {log} 1>&2"
