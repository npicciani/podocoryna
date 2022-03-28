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
		referencePeptides=expand("results/reference-isoseq/{transcriptome}.transdecoder_dir/longest_orfs.pep", transcriptome=config["reference"]["fileStem"]),
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

		ln -s {params.referencePeptides} {params.link}
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

# Provide species tree to orthofinder
# Include symbolic link to Podocoryna peptides
