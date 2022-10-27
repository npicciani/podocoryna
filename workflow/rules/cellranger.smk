rule cellranger_mkref_cogent:
    """
    Make reference using cell ranger and the original processed transcriptome file
    """
    input:
        fasta=expand("results/reference-cogent/{transcriptome}.fasta.selected.fixed", transcriptome=config["reference"]["fileStem"]),
        gtf=expand("results/reference-cogent/{transcriptome}.fasta.selected.fixed.eggnog.gtf", transcriptome=config["reference"]["fileStem"])
    output:
        directory("results/cellranger-cogent/reference")
    threads: 8
    shell:
        "cd results/cellranger-cogent && cellranger mkref --genome=reference --fasta=../../{input.fasta} --genes=../../{input.gtf}"

rule cellranger_count_cogent:
    input:
        "results/cellranger-cogent/reference"
    output:
        directory(expand("results/cellranger-cogent/counts/{species}", species=config["species"]))
    threads: 8
    params:
        species=config["species"]
    shell:
        """
        cd results/cellranger-cogent/counts
        cellranger count --id={params.species} \
                         --transcriptome=../reference \
                         --fastqs=../../../resources/rawdata \
                         --force-cells=8000 \
                         --localcores={threads} \
                         --localmem=64
        """
