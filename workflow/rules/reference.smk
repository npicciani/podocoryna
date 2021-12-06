rule generate_longest_ORFs:
    """
    Infer open reading frames in reference transcriptome.
    """
    input:
        transcriptomePath=expand("resources/{transcriptome}", transcriptome=config["reference"]["filename"]),
    output:
        expand("results/reference/{transcriptome}.transdecoder_dir/longest_orfs.pep", transcriptome=config["reference"]["filename"])
    params:
        referenceFilename=expand("{transcriptome}", transcriptome=config["reference"]["filename"])
    conda:
        "../../workflow/envs/transdecoder.yaml" #transdecoder v5.5.0
    shell:
        "TransDecoder.LongOrfs -t {input.transcriptomePath} --output_dir results/reference/{params.referenceFilename}.transdecoder_dir"

rule keep_longest_ORF_per_gene:
    """
    Keep longest open reading frames per gene (as per gene identifier type defined in the config yaml).
    Details on functions and arguments in the python script itself.
    """
    input:
        longestORFs=expand("results/reference/{transcriptome}.transdecoder_dir/longest_orfs.pep", transcriptome=config["reference"]["filename"]),
        transcriptomePath=expand("resources/{transcriptome}", transcriptome=config["reference"]["filename"]),
        script="workflow/scripts/keepLongestORFperGene.py"
    output:
        peptides=expand("results/reference/{transcriptome}_longestORFperGene.pep", transcriptome=config["reference"]["filename"]),
        nucleotides=expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"]["filename"])
    params:
        geneID_type=config["geneIDType"]
    shell:
        "python {input.script} -p {input.longestORFs} -t "
        "{input.transcriptomePath} -identifier {params.geneID_type} -o results/reference"

# rule make_GTF:
#     input:
#         nucleotides=expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"]["filename"]),
#         peptides=expand("results/reference/{transcriptome}_longestORFperGene.pep", transcriptome=config["reference"]["filename"]),
#         script="workflow/scripts/makeGTF_emapper.py"
#     output:
#         expand("results/reference/{transcriptome}_longestORFperGene.fasta.eggnog.gtf", transcriptome=config["reference"]["filename"]),
#         expand("results/reference/{transcriptome}_longestORFperGene.fasta.geneID_to_transcript.txt", transcriptome=config["reference"]["filename"])
#     threads: 15
#     params:
#         time="3:00:00",
#         mem="50GB",
#         geneID_type=config["geneIDType"]
#     shell:
#         "python {input.script} {input.nucleotides} {input.peptides} {params.geneID_type} results/reference"

rule get_mitochondrial_genes:
    input:
        nucleotides=expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"]["filename"]),
    output:
        directory("results/mitofinder")
    params:
        mem="50GB",
        mit_reference="resources/LN901210.1.gb"
    threads: 20
    shell:
        """
        mkdir {output} && cd {output}
        mitofinder -j podocoryna -a ../../{input.nucleotides} -r ../../{params.mit_reference} -o 4 -m {params.mem} -p {threads}
        """

rule get_annotated_transcripts:
    input:
        nucleotides=expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"]["filename"]),
        script="workflow/scripts/getSequences_onefasta.py"
    output:
        expand("results/reference/{transcriptome}_longestORFperGene.selected.fasta", transcriptome=config["reference"]["filename"])
    threads: 20
    params:
        transcriptlist="results/trinotate/transcripts_with_annotation_nodupes.txt"
    shell:
        """
        python {input.script} {params.transcriptlist} {input.nucleotides} results/reference
        """

rule select_from_gtf:
    input:
        gtf=expand("results/reference/{transcriptome}_longestORFperGene.fasta.eggnog.gtf", transcriptome=config["reference"]["filename"]),
        script="workflow/scripts/select_from_gtf.py"
    output:
        expand("results/reference/{transcriptome}_longestORFperGene.fasta.eggnog.selected.gtf", transcriptome=config["reference"]["filename"])
    params:
        transcriptlist="results/trinotate/transcripts_with_annotation.txt"
    shell:
        """
        python {input.script} {params.transcriptlist} {input.gtf} results/reference
        """
