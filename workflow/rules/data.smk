import pandas as pd

targets = (
    pd.read_csv(config["targets"], sep="\t", dtype={"species": str})
    .set_index("species", drop=False)
    .sort_index()
)

ensembl_targets = targets.loc[lambda targets: targets['source'] == "ensembl"]
ensemblgenomes_targets = targets.loc[lambda targets: targets['source'] == "ensemblgenomes"]
other_targets = targets.loc[lambda targets: targets['source'] == "other"]
other_gz_targets = targets.loc[lambda targets: targets['source'] == "other_gz"]

def get_sequence(wildcards, type):
    """
    Return path to download sequence file
    """
    species_units = targets.loc[wildcards.species]
    if type == "ensembl":
        if    species_units["source"] == "ensembl":
            return species_units["file"]
    if type == "ensemblgenomes":
        if    species_units["source"] == "ensemblgenomes":
            return species_units["file"]
    if type == "other":
        if species_units["source"] == "other":
            return species_units["file"]
    if type == "other_gz":
        if species_units["source"] == "other_gz":
            return species_units["file"]

rule get_ensembl:
    output:
        "resources/sequences/ensembl/{species}.pep.fasta.gz" # this renames the original file and keep only the species names as given in the species column of targets dataframe
    params:
        lambda wc: get_sequence(wc, type="ensembl")
    shell:
        "wget http://ftp.ensembl.org/{params} -O {output}"

rule get_ensemblgenomes:
    output:
        "resources/sequences/ensemblgenomes/{species}.pep.fasta.gz" # this renames the original file and keep only the species names as given in the species column of targets dataframe
    params:
        lambda wc: get_sequence(wc, type="ensemblgenomes")
    shell:
        "wget http://ftp.ensemblgenomes.org/{params} -O {output}"

rule wget_other:
    output:
        "resources/sequences/other/{species}.pep.fasta"
    params:
        lambda wc: get_sequence(wc, type="other")
    shell:
        """
        wget --no-check-certificate {params} -O {output}
        """

rule wget_other_gz:
    output:
        "resources/sequences/other_gz/{species}.pep.fasta.gz"
    params:
        lambda wc: get_sequence(wc, type="other_gz")
    shell:
        """
        wget --no-check-certificate {params} -O {output}
        """
