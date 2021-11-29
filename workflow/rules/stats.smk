rule busco_scores:
    input:
        expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"]["fileStem"])
    output:
        path=directory("results/busco")
    threads: 20
    conda:
        "../envs/busco.yaml" #busco=5.1.3
    params:
        download_path="results/busco/busco_downloads",
        mode="transcriptome",
        lineage="metazoa",
        filename=expand("{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"]["fileStem"])
    shell:
        "busco -i {input} -o {params.filename} --force --out_path {output.path} -l {params.lineage} -m {params.mode} --download_path {params.download_path} -c {threads}"
