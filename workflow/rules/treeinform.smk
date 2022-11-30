THRESHOLD_VALS = "0.5,1,2,3,4,5,8,10,15".split(",")

rule collapse_with_treeinform:
    input:
        gene_trees=expand("results/orthofinder/Results_{date}/Gene_Trees",date=ORTHODATE),
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
        gene_trees=expand("results/orthofinder/Results_{date}/Gene_Trees",date=ORTHODATE),
        script="workflow/scripts/plot_subtree_histogram.py",
        collapsed_proteins=expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.fasta", species=config["species"])
    output:
        "results/reference/treeinform/threshold_{threshold}/branch.length.hist_{threshold}.png"
    shell:
        """
        python {input.script} -gt {input.gene_trees} -t {wildcards.threshold} -b 20000 -o {output}
        """

rule match_transcripts:
    input:
        original_reference=expand("{transcriptome_path}", transcriptome_path=config["reference"]["path"]),
        collapsed_proteins=expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.fasta", species=config["species"]),
        script="workflow/scripts/pull_transcripts.py"
    output:
        transcripts=expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.fasta.transcripts.fasta", species=config["species"]),
        list_of_transcripts=expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.list.txt", species=config["species"])
    params:
        peptides=expand("results/reference/{transcriptome_stem}.transdecoder_dir/longest_orfs.pep", transcriptome_stem=config["reference"]["filestem"]),
    shell:
        """
        python {input.script} {input.original_reference} {input.collapsed_proteins}
        """

rule select_from_gtf:
    input:
        gtf=expand("{original_gtf}", original_gtf=config["reference"]["gtf"]),
        script="workflow/scripts/select_from_gtf.py",
        list_of_transcripts=expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.list.txt", species=config["species"])
    output:
        expand("results/reference/treeinform/threshold_{{threshold}}/{original_gtf_name}.selected.gtf", original_gtf_name=config["reference"]["gtfname"])
    shell:
        """
        python {input.script} {input.list_of_transcripts} {input.gtf} {output}
        """

# rule get_mitochondrial_genes:
#     input:
#         nucleotides=expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"]["filename"]),
#     output:
#         directory("results/mitofinder")
#     params:
#         mem="50GB",
#         mit_reference="resources/LN901210.1.gb",
#         mitofinder = expand("{mitofinderPath}", mitofinderPath=config["mitofinder"])
#     threads: 20
#     shell:
#         """
#         mkdir {output} && cd {output}
#         {params.mitofinder} -j podocoryna -a ../../{input.nucleotides} -r ../../{params.mit_reference} -o 4 -m {params.mem} -p {threads}
#         """
