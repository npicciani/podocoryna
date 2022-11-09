rule cellranger_mkref:
    """
    Make references using cell ranger and the thresholded transcriptome files
    """
    input:
        transcript=expand("results/reference/treeinform/threshold_{threshold}/{species}.collapsed.fasta.transcripts.fasta", threshold=["0.5","1","2","3","4","5","8","10","15"], species=config["species"]),
        gtf=expand("results/reference/treeinform/threshold_{threshold}/{original_gtf_name}.selected.gtf", threshold=["0.5","1","2","3","4","5","8","10","15"], original_gtf=config["reference"]["gtfname"])
    output:
        directory(expand("results/reference/treeinform/threshold_{threshold}/reference", threshold=["0.5","1","2","3","4","5","8","10","15"]))
    threads: 8
    params:
        outdir=expand("results/reference/treeinform/threshold_{threshold}", threshold=["0.5","1","2","3","4","5","8","10","15"])
    shell:
        "cd {params.outdir} && cellranger mkref --genome=reference --fasta=../../../../{input.transcript} --genes=../../../../{input.gtf}"

rule cellranger_count:
    input:
        directory(expand("results/reference/treeinform/threshold_{threshold}/reference", threshold=["0.5","1","2","3","4","5","8","10","15"]))
    output:
        directory(expand("results/reference/treeinform/threshold_{threshold}/{sample}", threshold=["0.5","1","2","3","4","5","8","10","15"], sample=["Podo1","Podo2","Podo3","Podo4","Podo5"]))
    threads: 8
    params:
        outdir=expand("results/reference/treeinform/threshold_{threshold}", threshold=["0.5","1","2","3","4","5","8","10","15"]),
        id=expand({sample}, sample=["Podo1","Podo2","Podo3","Podo4","Podo5"])
    shell:
        """
        cd {params.outdir}
        cellranger count --id={params.id} \
                         --transcriptome=./reference \
                         --fastqs=../../../../resources/rawdata/{params.id} \
                         --localcores={threads} \
                         --localmem=64
        """
