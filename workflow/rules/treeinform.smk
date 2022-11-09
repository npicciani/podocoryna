rule collapse_with_treeinform:
    input:
        gene_trees="results/orthofinder/Results/Gene_Trees",
        script="workflow/scripts/treeinform_collapse.py"
    output:
        expand("results/reference/treeinform/threshold_{threshold}", threshold=[0.5,1,2,3,4,5,8,10,15])
    threads: 20
    params:
        peptides=expand("results/reference/{transcriptome_stem}.transdecoder_dir/longest_orfs.pep", transcriptome_stem=config["reference"]["filestem"]),
        sp=config["species"]
    shell:
        """
        mkdir {output}
        python {input.script} -s {params.peptides} -gt {input.gene_trees} -t {threads} -sp {params.sp} -o {output}
        """

rule plot_subtree_histogram:
    input:
        gene_trees="results/orthofinder/Results/Gene_Trees",
        script="workflow/scripts/plot_subtree_histogram.py"
    output:
        expand("results/reference/treeinform/threshold_{threshold}", threshold=[0.5,1,2,3,4,5,8,10,15])
    params:
        sp=config["species"],
		thresh=expand("{thresholds}", threshold=[0.5,1,2,3,4,5,8,10,15])
    shell:
        """
        python {input.script} -gt {input.gene_trees} -t {params.thresh} -b 20000 -o {output}
        """

rule match_transcripts:
    input:
        original_reference=expand("{transcriptome_path}", transcriptome_path=config["reference"]["path"]),
        collapsed_proteins=expand("results/reference/treeinform/threshold_{threshold}/{species}.collapsed.fasta", threshold=[0.5,1,2,3,4,5,8,10,15], species=config["species"]),
        script="workflow/scripts/pull_transcripts.py"
    output:
        transcripts=expand("results/reference/treeinform/threshold_{threshold}/{species}.collapsed.fasta.transcripts.fasta", threshold=[0.5,1,2,3,4,5,8,10,15], species=config["species"]),
        list_of_transcripts=expand("results/reference/treeinform/threshold_{threshold}/{species}.collapsed.list.txt", threshold=[0.5,1,2,3,4,5,8,10,15], species=config["species"])
    threads: 20
    params:
        peptides=expand("results/reference/{transcriptome_stem}.transdecoder_dir/longest_orfs.pep", transcriptome_stem=config["reference"]["filestem"]),
        sp=config["species"]
    shell:
        """
        python {input.script} {input.original_reference} {input.collapsed_proteins}
        """

rule select_from_gtf:
    input:
        gtf=expand("{original_gtf}", original_gtf=config["reference"]["gtf"]),
        script="workflow/scripts/select_from_gtf.py"
    output:
        expand("results/reference/treeinform/threshold_{threshold}/{original_gtf_name}.selected.gtf", threshold=[0.5,1,2,3,4,5,8,10,15], original_gtf=config["reference"]["gtfname"])
    params:
        list_of_transcripts=expand("results/reference/treeinform/threshold_{threshold}/{species}.collapsed.list.txt", threshold=[0.5,1,2,3,4,5,8,10,15], species=config["species"])
    shell:
        """
        python {input.script} {params.list_of_transcripts} {input.gtf} {output}
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
