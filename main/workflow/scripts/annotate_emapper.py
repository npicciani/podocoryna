#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python 2.7
# -*- coding: utf-8 -*-
"""

Created by Jacob Musser
modified by Natasha Picciani
Jan 2023

Annotating proteins with eggNOG mapper, then summarizing emapper file, generating gene names

usage: annotate_emapper.py proteins.fasta outputDirectory

"""

import re
import subprocess
from pathlib import PurePosixPath
import sys
import argparse
from Bio import SeqIO

# Path to programs and files
emapper = "/home/picciani/local/bin/eggnog-mapper-2/emapper.py"  # path to eggNOG mapper version 2
python = "/home/picciani/local/bin/anaconda2/envs/py27/bin/python"  # path to python 2.7
gonames_file = "/home/picciani/local/datasets/go_terms_2019.txt"  # GO Terms IDs and their corresponding names from 2019 GO release


# User inputs
proteinFile = sys.argv[1]  # protein file
outDir = sys.argv[2]  # output directory
basename = PurePosixPath(proteinFile).stem

# Functionally annotating the protein sequences with eggNOG-mapper
subprocess.call(
    [
        python,
        emapper,
        "-i",
        proteinFile,
        "-m",
        "diamond",
        "-o",
        basename,
        "--cpu",
        "40",
        "--output_dir",
        outDir,
    ]
)

# # Parsing the emapper annotations file
# emapperFile = outDir + "/" + basename + ".emapper.annotations"

# emapper_proteinID_field = 0
# emapper_gene_name_field = 5
# emapper_gene_description_field = 21
# emapper_goterms_field = 6

# prot_name_dict = {}
# prot_description_dict = {}
# prot_list = []

# for line in open(emapperFile, 'r'):
# 	if line[0] == "#":
# 		continue
# 	else:
# 		line = line.strip('\n').split('\t')
# 		prot = line[emapper_proteinID_field]
# 		prot_list.append(prot)
# 		name = line[emapper_gene_name_field]
# 		prot_name_dict[prot] = name

# #names file
# names_outfile_name = outDir + "/" + basename + "_protein_names.txt"
# names_header = "proteinID\tshortname\tlongname"
# gene_shortname_dict = {}


# names_header = "geneID\tshortname\tlongname\tlongest_protein\tprotein_list"
# names_outlist = [names_header]
# for gene in gene_longest_protein_dict:
#     shortname = gene_shortname_dict[gene]
#     longname = gene_longname_dict[gene]
#     longest_protein = gene_longest_protein_dict[gene]
#     protein_list = geneID_proteinID_dict[gene]

#     names_outlist.append(gene + '\t' + shortname + '\t' + longname + '\t' + longest_protein + '\t' + ','.join(protein_list))

# names_outfile = open(names_outfile_name, 'w')
# for i in names_outlist:
#     names_outfile.write(i + '\n')
# names_outfile.close()


# # Creating final gene name and goterm dictionaries
# gene_shortname_dict = {}
# gene_longname_dict = {}
# gene_longest_protein_dict = {}
# gene_goterm_dict = {}


# for gene in geneID_proteinID_dict:
#     if len(geneID_proteinID_dict[gene]) > 0:
#         longest_protein = geneID_proteinID_dict[gene][0]
#         gene_longest_protein_dict[gene] = longest_protein

#         if longest_protein in prot_name_dict:
#         	if len(prot_name_dict[longest_protein]) > 0:
#         		gene_longname_dict[gene] = gene + " " + prot_name_dict[longest_protein] + " (" + prot_description_dict[longest_protein] +  ")"
#         		gene_shortname_dict[gene] = gene + " " + prot_name_dict[longest_protein]
#         	elif len(prot_description_dict[longest_protein]) >  0:
#         		gene_longname_dict[gene] = gene + " (" + prot_description_dict[longest_protein] +  ")"
#         		gene_shortname_dict[gene] = gene
#         	else:
#         		gene_longname_dict[gene] = gene
#         		gene_shortname_dict[gene] = gene
#         else:
#         	gene_longname_dict[gene] = gene
#         	gene_shortname_dict[gene] = gene

#         if longest_protein in prot_goterm_dict:
#         	if len(prot_goterm_dict[longest_protein][0]) > 0:
#         		gene_goterm_dict[gene] = prot_goterm_dict[longest_protein]
#         	else:
#         		gene_goterm_dict[gene]  = []
#         else:
#         	gene_goterm_dict[gene]  = []
#     else:
#     	gene_longest_protein_dict[gene] = ''
#     	gene_shortname_dict[gene] = gene
#     	gene_longname_dict[gene] = gene
