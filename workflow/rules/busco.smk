rule busco_scores:
    input:
        expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.fasta.transcripts.fasta", species=config["species"])
    output:
        directory("results/reference/treeinform/threshold_{threshold}/busco")
    threads: 20
    conda:
        "../envs/busco.yaml" #busco=5.1.3
    params:
        download_path="resources/busco_downloads",
        mode="transcriptome",
        lineage="metazoa",
        filename=expand("{species}.collapsed.fasta.transcripts.fasta", species=config["species"])
    shell:
        "busco -i {input} -o {params.filename} --force --out_path {output} -l {params.lineage} -m {params.mode} --download_path {params.download_path} -c {threads}"

# rule report_metrics:
#     input:
#         cellranger_mapping="results/reference/treeinform/threshold_{threshold}/cellranger/{sample}/outs/metrics_summary.csv",
#         busco_summary=expand("results/reference/treeinform/threshold_{{threshold}}/busco/{species}.collapsed.fasta.transcripts.fasta/short_summary.specific.metazoa_odb10.{species}.collapsed.fasta.transcripts.fasta.txt", species=config["species"])
#     output:
