rule collapse_with_treeinform:
    input:
        gene_trees=get_orthofinder_outdir(),
        script="workflow/scripts/treeinform_collapse.py"
    output:
        collapsed_proteins=expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.fasta", species=config["species"])
    params:
        peptides=expand("results/reference/{transcriptome_stem}.transdecoder_dir/longest_orfs.pep", transcriptome_stem=config["reference"]["filestem"]),
        sp=config["species"],
        outdir="results/reference/treeinform/threshold_{threshold}"
    shell:
        """
        mkdir {params.outdir}
        python {input.script} -s {params.peptides} -gt {input.gene_trees} -t {wildcards.threshold} -sp {params.sp} -o {params.outdir}
        """

rule plot_subtree_histogram:
    input:
        gene_trees=get_orthofinder_outdir(),
        script="workflow/scripts/plot_subtree_histogram.py",
        collapsed_proteins=expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.fasta", species=config["species"])
    output:
        image="results/reference/treeinform/threshold_{threshold}/branch.length.hist_{threshold}.png",
    params:
        outdir="results/reference/treeinform/threshold_{threshold}"
    shell:
        """
        python {input.script} -gt {input.gene_trees} -t {wildcards.threshold} -b 20000 -o {params.outdir}
        """

rule match_transcripts:
    input:
        original_reference=expand("{transcriptome_path}", transcriptome_path=config["reference"]["path"]),
        collapsed_proteins=expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.fasta", species=config["species"]),
        script="workflow/scripts/pull_transcripts.py"
    output:
        transcripts=expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.fasta.transcripts.fasta", species=config["species"]),
        list_of_transcripts=expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.fasta.transcripts.list.txt", species=config["species"])
    shell:
        """
        python {input.script} {input.original_reference} {input.collapsed_proteins}
        """

rule select_from_gtf:
    input:
        gtf=expand("{original_gtf}", original_gtf=config["reference"]["gtf"]),
        script="workflow/scripts/select_from_gtf.py",
        list_of_transcripts=expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.fasta.transcripts.list.txt", species=config["species"])
    output:
        expand("results/reference/treeinform/threshold_{{threshold}}/{original_gtf_name}.selected.gtf", original_gtf_name=config["reference"]["gtfname"])
    params:
        outdir="results/reference/treeinform/threshold_{threshold}"
    shell:
        """
        python {input.script} {input.list_of_transcripts} {input.gtf} {params.outdir}
        """
