rule makeblastdb:
    """
    Make blast database for reference from Matthew and Sally
    """
    input:
        expand("{proteome}", proteome=config["reference"]["proteome"])
    output:
        "results/blastx/Podocoryna_carnea.collapsed.fasta.psq"
    params:
        dbname="Podocoryna_carnea.collapsed.fasta"
    conda:
        "../envs/trinotate.yaml"
    shell:
        """
        cd results/blastx
        makeblastdb -in ../../{input} -out {params.dbname} -dbtype prot
        """
rule blastx:
    """
    Run blastx
    """
    input:
        "results/blastx/Podocoryna_carnea.collapsed.fasta.psq",
        fastafile = "resources/Podocoryna_lit_markers.fasta"
    output:
        "results/blastx/blastx.outfmt6"
    threads: 12
    conda:
        "../envs/trinotate.yaml"
    shell:
        """
        cd results/blastx

        # Blastx against reference db
        blastx \
        -query ../../{input.fastafile} \
        -db Podocoryna_carnea.collapsed.fasta \
        -num_threads 12 \
        -outfmt "6 qseqid sseqid qstart qend sstart send qstrand evalue" \
        -evalue 1e-8 > blastx.outfmt6
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
