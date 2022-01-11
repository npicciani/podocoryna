rule cellranger_mkref_original:
    """
    Make reference using cell ranger and the original processed transcriptome file
    """
    input:
        fasta=expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"]["filename"]),
        gtf=expand("results/reference/{transcriptome}_longestORFperGene.fasta.eggnog.gtf", transcriptome=config["reference"]["filename"])
    output:
        directory("results/cellranger/reference")
    threads: 8
    shell:
        "cd results/cellranger && cellranger mkref --genome=reference --fasta=../../{input.fasta} --genes=../../{input.gtf}"

rule cellranger_count_original:
    input:
        "results/cellranger/reference"
    output:
        directory(expand("results/cellranger/counts/{species}", species=config["species"]))
    threads: 8
    params:
        species=config["species"]
    shell:
        """
        cd results/cellranger/counts
        cellranger count --id={params.species} \
                         --transcriptome=../reference \
                         --fastqs=../../../resources/rawdata \
                         --expect-cells=10000 \
                         --localcores={threads} \
                         --localmem=64
        """

rule cellranger_mkref_selected:
    input:
        fasta=expand("results/reference/{transcriptome}_longestORFperGene.selected.fasta", transcriptome=config["reference"]["filename"]),
        gtf=expand("results/reference/{transcriptome}_longestORFperGene.fasta.eggnog.selected.gtf", transcriptome=config["reference"]["filename"])
    output:
        directory("results/cellranger/reference-selected")
    threads: 8
    shell:
        "cd results/cellranger && cellranger mkref --genome=reference-selected --fasta=../../{input.fasta} --genes=../../{input.gtf}"

rule cellranger_count_selected:
    input:
        "results/cellranger/reference-selected"
    output:
        directory(expand("results/cellranger/counts-selected/{species}", species=config["species"]))
    threads: 8
    params:
        species=config["species"]
    shell:
        """
        cd results/cellranger/counts-selected
        cellranger count --id={params.species} \
                         --transcriptome=../reference-selected \
                         --fastqs=../../../resources/rawdata \
                         --expect-cells=10000 \
                         --localcores={threads} \
                         --localmem=64
        """

rule cellranger_mkref_selected_annotated:
    input:
        fasta=expand("results/reference/{transcriptome}_longestORFperGene.allpositive.selected.fasta", transcriptome=config["reference"]["filename"]),
        gtf=expand("results/reference/{transcriptome}_longestORFperGene.fasta.eggnog.selected.gtf", transcriptome=config["reference"]["filename"])
    output:
        directory("results/cellranger/reference-selected-annotated")
    threads: 8
    shell:
        "cd results/cellranger && cellranger mkref --genome=reference-selected-annotated --fasta=../../{input.fasta} --genes=../../{input.gtf}"

rule cellranger_count_selected_annotated:
    input:
        "results/cellranger/reference-selected-annotated"
    output:
        directory(expand("results/cellranger/counts-selected-annotated/{species}", species=config["species"]))
    threads: 8
    params:
        species=config["species"]
    shell:
        """
        cd results/cellranger/counts-selected-annotated
        cellranger count --id={params.species} \
                         --transcriptome=../reference-selected-annotated \
                         --fastqs=../../../resources/rawdata \
                         --expect-cells=10000 \
                         --localcores={threads} \
                         --localmem=64
        """
