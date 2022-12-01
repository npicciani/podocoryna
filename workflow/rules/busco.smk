rule busco_scores:
    input:
        expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.fasta.transcripts.fasta", species=config["species"])
    output:
        outdir="results/reference/treeinform/threshold_{threshold}/busco"
    threads: 20
    conda:
        "../envs/busco.yaml" #busco=5.1.3
    params:
        download_path="resources/busco_downloads",
        mode="transcriptome",
        lineage="metazoa",
        filename=expand("{species}.collapsed.fasta.transcripts.fasta", species=config["species"])
    shell:
        "busco -i {input} -o {params.filename} --force --out_path {output.outdir} -l {params.lineage} -m {params.mode} --download_path {params.download_path} -c {threads}"
