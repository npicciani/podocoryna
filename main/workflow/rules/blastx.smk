rule makeblastdb:
    """
    Make blast database for reference from Matthew and Sally
    """
    input:
        sally = "results/reference-isoseq-original/Male_Med.ALL.Clustered.hq.CLEAN.cd_99.Min200.transdecoder_dir/longest_orfs.pep",
        matthew = "results/reference-ncbi/Podocoryna_carnea_transcriptome_2021_clean_NCBI_latest.fasta_longestORFperGene.pep"
    output:
        "results/blastx-isoseq/Male_Med.ALL.Clustered.hq.CLEAN.cd_99.Min200.longest_orfs.pep.psq",
        "results/blastx-ncbi/Podocoryna_carnea_transcriptome_2021_clean_NCBI_latest.fasta_longestORFperGene.pep.psq"
    params:
        dbname_matthew = "Podocoryna_carnea_transcriptome_2021_clean_NCBI_latest.fasta_longestORFperGene.pep",
        dbname_sally = "Male_Med.ALL.Clustered.hq.CLEAN.cd_99.Min200.longest_orfs.pep"
    conda:
        "../envs/trinotate.yaml"
    shell:
        """
        cd results/blastx-ncbi
        makeblastdb -in ../../{input.matthew} -out {params.dbname_matthew} -dbtype prot

        cd ../blastx-isoseq
        makeblastdb -in ../../{input.sally} -out {params.dbname_sally} -dbtype prot

        """
rule blastxIsoseq:
    """
    Run blastx
    """
    input:
        "results/blastx-isoseq/Male_Med.ALL.Clustered.hq.CLEAN.cd_99.Min200.longest_orfs.pep.psq",
        fastafile = "resources/Podocoryna_lit_markers.fasta"
    output:
        "results/blastx-isoseq/blastx.outfmt6"
    threads: 12
    conda:
        "../envs/trinotate.yaml"
    shell:
        """
        cd results/blastx-isoseq

        # Blastx against reference db

        blastx \
        -query ../../{input.fastafile} \
        -db Male_Med.ALL.Clustered.hq.CLEAN.cd_99.Min200.longest_orfs.pep \
        -num_threads 12 \
        -outfmt "6 qseqid sseqid qstart qend sstart send qstrand evalue" \
        -evalue 1e-5 > blastx.outfmt6
        """

rule blastx:
    """
    Run blastx
    """
    input:
        "results/blastx-ncbi/Podocoryna_carnea_transcriptome_2021_clean_NCBI_latest.fasta_longestORFperGene.pep.psq",
        fastafile = "resources/Podocoryna_lit_markers.fasta"
    output:
        "results/blastx-ncbi/blastx.outfmt6"
    threads: 12
    conda:
        "../envs/trinotate.yaml"
    shell:
        """
        cd results/blastx-ncbi

        # Blastx against reference db

        blastx \
        -query ../../{input.fastafile} \
        -db Podocoryna_carnea_transcriptome_2021_clean_NCBI_latest.fasta_longestORFperGene.pep \
        -num_threads 12 \
        -outfmt "6 qseqid sseqid qstart qend sstart send qstrand evalue" \
        -evalue 1e-5 > blastx.outfmt6
        """

# rule makeblastdb:
#     """
#     Make Uniprot (trEmbl + swissprot) blast database
#     """
#     output:
#         "results/blastx/uniprot-taxonomy%3A6073.fasta.psq"
#     conda:
#         "../envs/trinotate.yaml"
#     shell:
#         """
#         cd results/blastx
#         makeblastdb -in uniprot-taxonomy%3A6073.fasta -dbtype prot
#         """
# rule blastx:
#     """
#     Run blastx
#     """
#     input:
#         "results/blastx/uniprot-taxonomy%3A6073.fasta",
#         "results/blastx/uniprot-taxonomy%3A6073.fasta.psq",
#         fastafile = expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"]["filename"]),
#     output:
#         "results/blastx/blastx.outfmt6"
#     threads: 12
#     conda:
#         "../envs/trinotate.yaml"
#     shell:
#         """
#         cd results/blastx

#         # Blastx against uniprot

#         blastx \
#         -query ../../{input.fastafile} \
#         -db uniprot-taxonomy%3A6073.fasta \
#         -num_threads 12 \
#         -max_target_seqs 1 \
#         -outfmt "6 qseqid sseqid qstart qend sstart send qstrand evalue" \
#         -evalue 1e-6 > blastx.outfmt6
#         """
