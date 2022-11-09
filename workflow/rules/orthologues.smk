rule generate_longest_ORFs:
    """
    Generate open reading frames from reference transcriptome.
    """
    input:
        transcriptome_path=expand("{transcriptome_path}", transcriptome_path=config["reference"]["path"])
    output:
        expand("results/reference/{transcriptome_stem}.transdecoder_dir/longest_orfs.pep", transcriptome_stem=config["reference"]["filestem"])
    params:
        reference=expand("{transcriptome_stem}", transcriptome_stem=config["reference"]["filestem"])
    conda:
        "../../workflow/envs/transdecoder.yaml" #transdecoder v5.5.0
    shell:
        "TransDecoder.LongOrfs -t {input.transcriptome_path} --output_dir results/reference/{params.reference}.transdecoder_dir"

rule gunzip:
	input:
		expand("resources/sequences/ensembl/{species}.pep.fasta.gz",species=ensembl_targets.loc[:,"species"]),
		expand("resources/sequences/ensemblgenomes/{species}.pep.fasta.gz",species=ensemblgenomes_targets.loc[:,"species"]),
		expand("resources/sequences/gdrive/{species}.pep.fasta",species=gdrive_targets.loc[:,"species"]),
		expand("resources/sequences/other/{species}.pep.fasta",species=other_targets.loc[:,"species"]),
		expand("resources/sequences/other_gz/{species}.pep.fasta.gz",species=other_gz_targets.loc[:,"species"])
	output:
		expand("resources/sequences/{species}.pep.fasta", species=targets.index)
	params:
		reference_peptides=expand("results/reference/{transcriptome_stem}.transdecoder_dir/longest_orfs.pep", transcriptome_stem=config["reference"]["filestem"]),
		link=expand("resources/sequences/{species}.pep.fasta", species=config["reference"]["species"])
	shell:
		"""
		dir="resources/sequences"
		gzfiles=`ls $dir/*/*.gz`
		for file in $gzfiles; do
			gunzip $file
		done
		fastafiles=`ls $dir/*/*.fasta`
		for file in $fastafiles; do
			mv $file $dir
		done
		subdirs=`ls -d $dir/*/`
		rm -R $subdirs

		ln -s {params.reference_peptides} {params.link}
		"""

rule orthofinder:
    input:
        directory("resources/sequences")
    output:
        directory("results/orthofinder")
    conda:
        "../envs/orthofinder.yaml"
    log:
        "logs/orthofinder/orthofinder.log"
    threads: 20
    shell:
        "orthofinder -f {input} -t {threads} -o {output}"
