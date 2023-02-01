from datetime import date
def get_orthofinder_outdir():
    """
    Generate path to orthofinder gene trees folder with current date as written by orthofinder
    """
    today = date.today()
    monthDay = today.strftime("%b%d")
    outdir= f"results/orthofinder/Results_{monthDay}/Gene_Trees"
    return outdir

rule gunzip:
    """
    Decompress (if gzipped) protein files downloaded from public databases.
    """
    input:
        expand("resources/sequences/ensembl/{species}.pep.fasta.gz",species=ensembl_targets.loc[:,"species"]),
        expand("resources/sequences/ensemblgenomes/{species}.pep.fasta.gz",species=ensemblgenomes_targets.loc[:,"species"]),
        expand("resources/sequences/other/{species}.pep.fasta",species=other_targets.loc[:,"species"]),
        expand("resources/sequences/other_gz/{species}.pep.fasta.gz",species=other_gz_targets.loc[:,"species"]),
        expand("resources/sequences/gdrive/{species}.pep.fasta",species=gdrive_targets.loc[:,"species"]),
        reference_peptides="../ref_optimization/results/reference/Male_Med.ALL.Clustered.hq.CLEAN.cd_99.Min200.transdecoder_dir/longest_orfs.pep"
    output:
        expand("resources/sequences/{species}.pep.fasta", species=targets.index)
    params:
        copyfile=expand("resources/sequences/{species}.pep.fasta", species=config["species"]),
        nanomia_proteins="/home/nnp9/scratch60/nanomia/PO_YAL3284_Nanomia_bijuga.protein.fasta"
    shell:
        """
        dir="resources/sequences"
        gzfiles=`find $dir/* -name '*.gz'`
        for file in $gzfiles; do
            gunzip $file
        done
        fastafiles=`find $dir/* -name '*.fasta'`
        for file in $fastafiles; do
            mv $file $dir
        done
        subdirs=`ls -d $dir/*/`
        rm -R $subdirs

        cp {input.reference_peptides} {params.copyfile}
        cp {params.nanomia_proteins} resources/sequences/Nanomia_bijuga.pep.fasta
        """

rule emapper_annotate:
    """
    Generate functional annotations using eggnog mapper.
    """
    input:
        expand("resources/sequences/{species}.pep.fasta", species=targets.index)
    output:
        "results/annotations/all.pep.fasta.emapper.annotations"
    params:
        outdir="results/annotations",
        script="workflow/scripts/annotate_emapper.py",
        peptides="resources/sequences/all.pep.fasta"
    conda:
        "../../workflow/envs/emapper.yaml" #eggnog-mapper=2.0.6, python=3.7.9
    shell:
        """
        cat resources/sequences/*.fasta > {params.peptides}
        python {params.script} {params.peptides} {params.outdir}
        """

rule orthofinder:
    """
    Infer gene trees from set of protein sequences downloaded from public databases.
    """
    input:
        expand("resources/sequences/{species}.pep.fasta", species=targets.index)
    output:
        directory(get_orthofinder_outdir())
    conda:
        "../../workflow/envs/orthofinder.yaml"
    threads: 20
    shell:
        """
        rm -rf results/orthofinder
        orthofinder -t {threads} -f resources/sequences -o results/orthofinder
        """
