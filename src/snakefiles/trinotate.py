rule trinotate_split_assembly:
    """
    Split the headers from assembly into multiple files
    """
    input:
        fai = raw + "assembly.fasta.fai"
    output:
        expand(
            trinotate + "chunks/assembly_{chunk_id}.tsv",
            chunk_id=[
                '{0:05d}'.format(x)
                for x in range(0, config["number_of_chunks"]["trinotate"])
            ]
        )
    params:
        number_of_chunks = config["number_of_chunks"]["trinotate"]
    log:
        trinotate + "split_assembly.log"
    benchmark:
        trinotate + "split_assembly.json"
    conda:
        "trinotate.yml"
    shell:
        """
        split \
            --number l/{params.number_of_chunks} \
            --numeric-suffixes \
            --suffix-length 5 \
            --additional-suffix .tsv \
            {input.fai} \
            {trinotate}/chunks/assembly_ \
        2> {log}
        """


rule trinotate_blastx_chunk:
    """
    Run blastx of each chunk
    """
    input:
        fa = raw + "assembly.fasta",
        fai = raw + "assembly.fasta.fai",
        chunk = trinotate + "chunks/assembly_{chunk_id}.tsv",
        db = db + "uniprot_sprot"
    output:
        tsv = trinotate + "blastx/assembly_{chunk_id}.tsv"
    log:
        trinotate + "blastx/assembly_{chunk_id}.log"
    benchmark:
        trinotate + "blastx/assembly_{chunk_id}.json"
    conda:
        "trinotate.yml"
    shell:
        """
        cut -f 1 {input.chunk} \
        | xargs samtools faidx {input.fa} \
        | blastx \
            -db {input.db} \
            -max_target_seqs 1 \
            -outfmt 6 \
            -evalue 1e-5 \
            -out {output.tsv} \
        2> {log} 1>&2
        """


rule trinotate_blastx_merge:
    """
    Merge results from the different blastxs
    """
    input:
        expand(
            trinotate + "blastx/assembly_{chunk_id}.tsv",
            chunk_id=[
                '{0:05d}'.format(x)
                for x in range(0, config["number_of_chunks"]["trinotate"])
            ]
        )
    output:
        tsv = trinotate + "blastx.tsv"
    log:
        trinotate + "blastx_merge.log"
    benchmark:
        trinotate + "blastx_merge.json"
    conda:
        "trinotate.yml"
    shell:
        "cat {input} > {output} 2> {log}"


rule trinotate_split_proteome:
    """
    Split the headers from transcriptome.pep into multiple files
    """
    input:
        fai = transdecoder + "transdecoder.pep.fai"
    output:
        expand(
            trinotate + "chunks/proteome_{chunk_id}.tsv",
            chunk_id=[
                '{0:05d}'.format(x)
                for x in range(0, config["number_of_chunks"]["trinotate"])
            ]
        )
    params:
        number_of_chunks = config["number_of_chunks"]["trinotate"]
    log:
        trinotate + "split_proteome.log"
    benchmark:
        trinotate + "split_proteome.json"
    conda:
        "trinotate.yml"
    shell:
        """
        split \
            --number l/{params.number_of_chunks} \
            --numeric-suffixes \
            --suffix-length 5 \
            --additional-suffix .tsv \
            {input.fai} \
            {trinotate}/chunks/proteome_ \
        2> {log}
        """


rule trinotate_blastp_chunk:
    """
    Run blastp of each chunk
    """
    input:
        pep = transdecoder + "transdecoder.pep",
        fai = transdecoder + "transdecoder.pep.fai",
        chunk = trinotate + "chunks/proteome_{chunk_id}.tsv",
        db = db + "uniprot_sprot"
    output:
        tsv = trinotate + "blastp/proteome_{chunk_id}.tsv"
    log:
        trinotate + "blastp/proteome_{chunk_id}.log"
    benchmark:
        trinotate + "blastp/proteome_{chunk_id}.json"
    conda:
        "trinotate.yml"
    shell:
        """
        cut -f 1 {input.chunk} \
        | xargs samtools faidx {input.pep} \
        | blastp \
            -db {input.db} \
            -max_target_seqs 1 \
            -outfmt 6 \
            -evalue 1e-5 \
            -out {output.tsv} \
        2> {log} 1>&2
        """


rule trinotate_blastp_merge:
    """
    Merge results from the different blastps
    """
    input:
        expand(
            trinotate + "blastp/proteome_{chunk_id}.tsv",
            chunk_id=[
                '{0:05d}'.format(x)
                for x in range(0, config["number_of_chunks"]["trinotate"])
            ]
        )
    output:
        tsv = trinotate + "blastp.tsv"
    log:
        trinotate + "blastp_merge.log"
    benchmark:
        trinotate + "blastp_merge.json"
    conda:
        "trinotate.yml"
    shell:
        "cat {input} > {output} 2> {log}"


rule trinotate_hmmscan_chunk:
    """
    hmmscan over one chunk
    """
    input:
        pep = transdecoder + "transdecoder.pep",
        fai = transdecoder + "transdecoder.pep.fai",
        chunk = trinotate + "chunks/proteome_{chunk_id}.tsv",
        hmm = db + "Pfam-A.hmm"
    output:
        tsv = trinotate + "hmmscan/proteome_{chunk_id}.tsv"
    log:
        trinotate + "hmmscan/proteome_{chunk_id}.log"
    benchmark:
        trinotate + "hmmscan/proteome_{chunk_id}.json"
    conda:
        "trinotate.yml"
    shell:
        """
        cut -f 1 {input.chunk} \
        | xargs samtools faidx {input.pep} \
        | hmmscan \
            --domtblout {output.tsv} \
            {input.hmm} \
            - \
        2> {log} 1>&2
        """


rule trinotate_hmmscan_merge:
    """
    Merge hmmscan results into one file
    """
    input:
        expand(
            trinotate + "hmmscan/proteome_{chunk_id}.tsv",
            chunk_id=[
                '{0:05d}'.format(x)
                for x in range(0, config["number_of_chunks"]["trinotate"])
            ]
        )
    output:
        tsv = trinotate + "hmmscan.tsv"
    log:
        trinotate + "hmmscan_merge.log"
    benchmark:
        trinotate + "hmmscan_merge.json"
    conda:
        "trinotate.yml"
    shell:
        "cat {input} > {output} 2> {log}"


# rule trinotate_signalp:
#     """
#     Use SignalP to predict signal peptides
#     """
#     input:
#         pep = transdecoder + "transdecoder.pep"
#     output:
#         tsv = touch(trinotate + "signalp.tsv")
#     log:
#         trinotate + "signalp.log"
#     benchmark:
#         trinotate + "signalp.json"
#     conda:
#         "trinotate.yml"
#     shell:
#         "./src/signalp-4.1/signalp "
#             "-f short "
#             "-n {output.tsv} "
#             "{input.pep} "
#         "2> {log} 1>&2 "


# rule trinotate_tmhmm:
#     """
#     Predict transmembrane regions
#     """
#     input:
#         pep = transdecoder + "transdecoder.pep"
#     output:
#         tsv = trinotate + "tmhmm.tsv"
#     log:
#         trinotate + "tmhmm.log"
#     benchmark:
#         trinotate + "tmhmm.json"
#     conda:
#         "trinotate.yml"
#     shell:
#         "./src/tmhmm-2.0c/bin/tmhmm "
#             "--short "
#         "< {input.pep} "
#         "> {output.tsv} "
#         "2> {log} && "
#         "rm -rf TMHMM_*"


# rule trinotate_rnammer:
#     """
#     Identify rRNAs
#     """
#     input:
#         assembly = raw + "assembly.fasta",
#     output:
#         gff = trinotate + "rnammer.gff"
#         "assembly.fasta.rnammer.gff"
#     params:
#         rnammer_path = config["trinotate"]["rnammer"]["rnammer_path"],
#         org_type = config["trinotate"]["rnammer"]["org_type"],
#         gff_tmp = "assembly.fasta.rnammer.gff"
#     log:
#         trinotate + "rnammer.log"
#     benchmark:
#         trinotate + "rnammer.json"
#     conda:
#         "trinotate.yml"
#     shell:
#         "./src/Trinotate-3.0.2/util/rnammer_support/RnammerTranscriptome.pl "
#             "--transcriptome {input.assembly} "
#             "--path_to_rnammer  ./src/rnammer-1.2/rnammer "
#             "--org_type {params.org_type} "
#         "2> {log} 1>&2 ; "
#         "mv {params.gff_tmp} {output.gff} ; "
#         "rm tmp.superscaff.rnammer.gff "
#             "transcriptSuperScaffold.bed "
#             "transcriptSuperScaffold.fasta"


rule trinotate_create:
    """
    Build sqlite database
    """
    output:
        sqlite = trinotate + "trinotate.sqlite",
        is_created = touch(trinotate + "create.txt")
    log:
        trinotate + "create.log"
    benchmark:
        trinotate + "create.json"
    conda:
        "trinotate.yml"
    shell:
        """
        EMBL_dat_to_Trinotate_sqlite_resourceDB.pl \
            --sqlite {output.sqlite} \
            --create \
        2> {log} 1>&2
        """


rule trinotate_init:
    """
    Initialize db with genes, transcripts and proteins
    """
    input:
        sqlite = ancient(trinotate + "trinotate.sqlite"),
        is_created = trinotate + "create.txt",
        g2t = raw + "gene_to_trans_map.tsv",
        assembly = raw + "assembly.fasta",
        proteome = transdecoder + "transdecoder.pep"
    output:
        touch(trinotate + "init.txt")
    log:
        trinotate + "init.log"
    benchmark:
        trinotate + "init.json"
    conda:
        "trinotate.yml"
    shell:
        """
        Trinotate {input.sqlite} init \
            --gene_trans_map {input.g2t} \
            --transcript_fasta {input.assembly} \
            --transdecoder_pep {input.proteome} \
        2> {log} 1>&2
        """


rule trinotate_fill:
    input:
        sqlite = ancient(trinotate + "trinotate.sqlite"),
        is_init = trinotate + "init.txt",
        eggnog = db + "NOG.annotations.tsv.bulk_load",
        go = db + "go-basic.obo.tab",
        uniprot_index = db + "trinotate.UniprotIndex",
        taxonomy_index = db + "trinotate.TaxonomyIndex",
        pfam = db + "Pfam-A.hmm.gz.pfam_sqlite_bulk_load"
    output:
        is_filled = touch(trinotate + "fill.txt")
    log:
        db + "fill.log"
    benchmark:
        db + "fill.json"
    conda:
        "trinotate.yml"
    shell:
        """
        EMBL_dat_to_Trinotate_sqlite_resourceDB.pl \
            --sqlite {input.sqlite} \
            --eggnog {input.eggnog} \
            --go_obo_tab {input.go} \
            --uniprot_index {input.uniprot_index} \
            --taxonomy_index {input.taxonomy_index} \
            --pfam {input.pfam} \
        2> {log} 1>&2
        """


rule trinotate_load:
    input:
        sqlite = ancient(trinotate + "trinotate.sqlite"),
        is_filled = trinotate + "fill.txt",
        blastx = trinotate + "blastx.tsv",
        blastp = trinotate + "blastp.tsv",
        pfam = trinotate + "hmmscan.tsv",
        # signalp = trinotate + "signalp.tsv",
        # tmhmm = trinotate + "tmhmm.tsv",
        # rnammer = trinotate + "rnammer.gff",
    output:
        touch(trinotate + "load.txt")
    log:
        trinotate + "load.log"
    benchmark:
        trinotate + "load.json"
    conda:
        "trinotate.yml"
    shell:
        """
        Trinotate \
            {input.sqlite} \
            LOAD_swissprot_blastp \
            {input.blastp} \
        2> {log} 1>&2

        Trinotate \
            {input.sqlite} \
            LOAD_pfam \
            {input.pfam} \
        2>> {log} 1>&2

        Trinotate \
            {input.sqlite} \
            LOAD_swissprot_blastx \
            {input.blastx} \
        2>> {log} 1>&2
        """
        # """
        # Trinotate \
        #     {input.sqlite} \
        #     LOAD_signalp \
        #     {input.signalp} \
        # 2>> {log} 1>&2
        #
        # Trinotate \
        #     {input.sqlite} \
        #     LOAD_tmhmm \
        #     {input.tmhmm} \
        # 2>> {log} 1>&2
        #
        # Trinotate \
        #   {input.sqlite} \
        #   LOAD_rnammer \
        #   {input.rnammer} \
        # 2>> {log} 1>&2
        # """


rule trinotate_report:
    input:
        sqlite = ancient(trinotate + "trinotate.sqlite"),
        is_loaded = trinotate + "load.txt"
    output:
        trinotate + "trinotate.tsv"
    params:
        evalue = config["trinotate"]["evalue"],
        pfam_cutoff = config["trinotate"]["pfam_cutoff"]
    log:
        trinotate + "report.log"
    benchmark:
        trinotate + "report.json"
    conda:
        "trinotate.yml"
    shell:
        """
        Trinotate {input.sqlite} report \
            -E {params.evalue} \
            --pfam_cutoff {params.pfam_cutoff} \
        > {output} \
        2> {log}
        """
