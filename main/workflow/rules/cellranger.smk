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
        transcript=config["reference"]["transcriptome_with_mito"],
        gtf=expand("resources/reference/{gtf}",gtf=config["reference"]["gtfname_with_mito"])
    output:
        directory("results/cellranger/reference")
    threads: 8
    params:
        outdir="results/cellranger"
    shell:
        """
        cd {params.outdir}
        cellranger mkref --genome=reference \
                        --fasta=../../{input.transcript} \
                        --genes=../../{input.gtf}
        """

rule count_cellranger:
    """
    Run cellranger count to align single cell reads to reference transcriptome.
    """

    input:
        input_dir="results/cellranger/reference"
    output:
        outfile="results/cellranger/{sample}/outs/molecule_info.h5",
    threads: 8
    params:
        outdir="results/cellranger",
        reference_dir="./reference",
        fastqs_dir="../../../ref_optimization/resources/rawdata/{sample}"
    shell:
        """
        rm -rf results/cellranger/{wildcards.sample}
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
        info="results/cellranger/{sample}/outs/molecule_info.h5"
    output:
        aggregate_init_stamp="results/cellranger/{sample}/{sample}.aggregate_init.stamp",
    params:
        relative_path="../{sample}/outs/molecule_info.h5",
        csv = "results/cellranger/aggregate/samples.csv",
        batch = get_batch
    run:
        import os
        if os.path.exists(params.csv):
            with open(params.csv,"a") as out1:
                out1.write(f"{(wildcards.sample)},{(params.relative_path)},{(params.batch)}\n")
        if not os.path.exists(params.csv):
            shell("mkdir results/cellranger/aggregate")
            with open(params.csv, "w") as out2:
                out2.write("sample_id,molecule_h5,batch"+"\n")
                out2.write(f"{(wildcards.sample)},{(params.relative_path)},{(params.batch)}\n")
        shell("touch {output.aggregate_init_stamp}")

rule aggregate_cellranger:
    """
    Aggregate count data from multiple samples using cellranger aggregate.
    """
    input:
        expand("results/cellranger/{sample}/{sample}.aggregate_init.stamp", sample=SAMPLE_IDS)
    output:
        aggregate="results/cellranger/aggregate/aggregated/outs/web_summary.html"
    params:
        outdir = "results/cellranger/aggregate"
    shell:
        """
        rm -rf results/cellranger/aggregate/aggregated
        cd {params.outdir}
        cellranger aggr --id=aggregated \
                        --csv=samples.csv
        """
