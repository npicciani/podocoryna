def get_batch(wildcards):
    """
    Return string with batch information from each sample
    """
    return config["batch"][wildcards.sample]

rule mkref_cellranger:
    """
    Make references using cellranger and the thresholded transcriptome files
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

rule count_cellranger:
    """
    Run cellranger count to align single cell reads to reference transcriptome.
    """

    input:
        input_dir="results/reference/treeinform/threshold_{threshold}/cellranger/reference"
    output:
        outfile="results/reference/treeinform/threshold_{threshold}/cellranger/{sample}/outs/molecule_info.h5",
    threads: 8
    params:
        outdir="results/reference/treeinform/threshold_{threshold}/cellranger",
        reference_dir="./reference",
        fastqs_dir="../../../../../resources/rawdata/{sample}"
    shell:
        """
        rm -rf results/reference/treeinform/threshold_{wildcards.threshold}/cellranger/{wildcards.sample}
        cd {params.outdir}
        cellranger count --id={wildcards.sample} \
                         --transcriptome={params.reference_dir} \
                         --fastqs={params.fastqs_dir} \
                         --localcores={threads} \
                         --localmem=64
        """
rule aggregate_init_cellranger:
    """
    Start and populate aggregation CSV file with sample information that
    is input for cellranger aggregate. Generate a stamp for each sample once
    its info is added to CSV file.
    """
    input:
        info="results/reference/treeinform/threshold_{threshold}/cellranger/{sample}/outs/molecule_info.h5"
    output:
        aggregate_init_stamp="results/reference/treeinform/threshold_{threshold}/cellranger/{sample}/{sample}.aggregate_init.stamp",
    params:
        relative_path="../../cellranger/{sample}/outs/molecule_info.h5",
        csv = "results/reference/treeinform/threshold_{threshold}/cellranger/aggregate/samples.csv",
        batch = get_batch
    run:
        import os
        if os.path.exists(params.csv):
            with open(params.csv,"a") as out1:
                out1.write(f"{(wildcards.sample)},{(params.relative_path)},{(params.batch)}\n")
        if not os.path.exists(params.csv):
            shell("mkdir results/reference/treeinform/threshold_{wildcards.threshold}/cellranger/aggregate")
            with open(params.csv, "w") as out2:
                out2.write("sample_id,molecule_h5,batch"+"\n")
                out2.write(f"{(wildcards.sample)},{(params.relative_path)},{(params.batch)}\n")
        shell("touch {output.aggregate_init_stamp}")

rule aggregate_cellranger:
    """
    Aggregate count data from multiple samples using cellranger aggregate.
    """

    input:
        aggregate_init_stamp="results/reference/treeinform/threshold_{threshold}/cellranger/{sample}/{sample}.aggregate_init.stamp"
    output:
        aggregate_stamp="results/reference/treeinform/threshold_{threshold}/cellranger/{sample}/{sample}.aggregate.stamp"
    params:
        outdir = "results/reference/treeinform/threshold_{threshold}/cellranger/aggregate"
    shell:
        """
        cd {params.outdir}
        cellranger aggr --id=aggregated \
                        --csv=samples.csv
        touch ../../../../../../{output.aggregate_stamp}
        """
