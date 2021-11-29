rule cellranger_mkref:
    input:
        fasta=expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"]["filename"]),
        gtf=expand("results/reference/{transcriptome}_longestORFperGene.fasta.eggnog.gtf", transcriptome=config["reference"]["filename"]),
    output:
        directory("results/cellranger/reference")
    threads: 8
    shell:
        "cd results/cellranger && cellranger mkref --genome=reference --fasta=../../{input.fasta} --genes=../../{input.gtf}"

rule cellranger_count:
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