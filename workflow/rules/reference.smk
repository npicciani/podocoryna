rule generate_longest_ORFs:
    """
    Infer open reading frames in reference transcriptome.
    """
    input:
        transcriptomePath=expand("results/reference-cogent/{transcriptome}.fixed", transcriptome=config["reference"]["filename"])
    output:
        expand("results/reference-cogent/{transcriptome}.fixed.transdecoder_dir/longest_orfs.pep", transcriptome=config["reference"]["filename"])
    params:
        referenceFilename=expand("{transcriptome}", transcriptome=config["reference"]["filename"])
    conda:
        "../../workflow/envs/transdecoder.yaml" #transdecoder v5.5.0
    shell:
        "TransDecoder.LongOrfs -t {input.transcriptomePath} --output_dir results/reference-cogent/{params.referenceFilename}.fixed.transdecoder_dir"


# rule get_mitochondrial_genes:
#     input:
#         nucleotides=expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"]["filename"]),
#     output:
#         directory("results/mitofinder")
#     params:
#         mem="50GB",
#         mit_reference="resources/LN901210.1.gb",
#         mitofinder = expand("{mitofinderPath}", mitofinderPath=config["mitofinder"])
#     threads: 20
#     shell:
#         """
#         mkdir {output} && cd {output}
#         {params.mitofinder} -j podocoryna -a ../../{input.nucleotides} -r ../../{params.mit_reference} -o 4 -m {params.mem} -p {threads}
#         """

# rule get_annotated_transcripts:
#     input:
#         nucleotides=expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"]["filename"]),
#         script="workflow/scripts/getSequences_onefasta.py"
#     output:
#         expand("results/reference/{transcriptome}_longestORFperGene.selected.fasta", transcriptome=config["reference"]["filename"])
#     threads: 20
#     params:
#         transcriptlist="results/trinotate/transcripts_with_annotation_nodupes.txt"
#     shell:
#         """
#         python {input.script} {params.transcriptlist} {input.nucleotides} results/reference
#         """

# rule select_from_gtf:
#     input:
#         gtf=expand("results/reference/{transcriptome}_longestORFperGene.fasta.eggnog.gtf", transcriptome=config["reference"]["filename"]),
#         script="workflow/scripts/select_from_gtf.py"
#     output:
#         expand("results/reference/{transcriptome}_longestORFperGene.fasta.eggnog.selected.gtf", transcriptome=config["reference"]["filename"])
#     params:
#         transcriptlist="results/trinotate/transcripts_with_annotation.txt"
#     shell:
#         """
#         python {input.script} {params.transcriptlist} {input.gtf} results/reference
#         """
