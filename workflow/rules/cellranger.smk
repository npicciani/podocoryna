rule cellranger_mkref:
    """
    Make references using cell ranger and the thresholded transcriptome files
    """
    input:
        transcript=expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.fasta.transcripts.fasta", species=config["species"]),
        gtf=expand("results/reference/treeinform/threshold_{{threshold}}/{original_gtf_name}.selected.gtf", original_gtf_name=config["reference"]["gtfname"])
    output:
        directory("results/reference/treeinform/threshold_{threshold}/cellranger/reference")
    threads: 8
    params:
        outdir="results/reference/treeinform/threshold_{threshold}/cellranger"
    shell:
        """
        cd {params.outdir}
        cellranger mkref --genome=reference \
                        --fasta=../../../../../{input.transcript} \
                        --genes=../../../../../{input.gtf}
        """

rule cellranger_count:
    input:
        input_dir="results/reference/treeinform/threshold_{threshold}/cellranger/reference"
    output:
        outfile="results/reference/treeinform/threshold_{threshold}/cellranger/{sample}/outs/web_summary.html"
    threads: 8
    params:
        outdir="results/reference/treeinform/threshold_{threshold}/cellranger",
        id="{sample}",
        reference_dir="./reference",
        fastqs_dir="../../../../../resources/rawdata/{sample}"
    shell:
        """
        rm -rf results/reference/treeinform/threshold_{wildcards.threshold}/cellranger/{wildcards.sample}
        cd {params.outdir}
        cellranger count --id={params.id} \
                         --transcriptome={params.reference_dir} \
                         --fastqs={params.fastqs_dir} \
                         --localcores={threads} \
                         --localmem=64
        """
