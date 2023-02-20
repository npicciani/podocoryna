# -*- coding: utf-8 -*-
# Written by Natasha Picciani
# Oct 31, 2022
# usage: identify_transcripts.py fastafile

import re
import sys

annotations_master = sys.argv[1]
gene_trees_master = sys.argv[2]
outdir = sys.argv[3]
gene_trees_appended = f"{outdir}/gene_trees.master.annotated.txt"

gene_names_dict = {}
# GOs_dict={}
with open(annotations_master, "r") as infile:
    for line in infile:
        if line[0] == "#":
            continue
        else:
            line = line.strip("\n")
            columns = line.split("\t")
            proteinID = columns[0]
            proteinID_gene_name = f"{columns[0]}_{columns[8]}"
            gene_name = columns[8]  # Preferred_name
            gene_description = columns[7]  # best_og_desc
            # goterms_field = columns[9]  # GOs
            gene_names_dict[proteinID] = proteinID_gene_name
            # GOs_dict[proteinID] = goterms_field
    # outfile.write(proteinID_gene_name+"\t"+proteinID+"\t"+gene_name+"\t"+gene_description+"\t"+goterms_field+"\n")

with open(gene_trees_master) as trees:
    lines = trees.readlines()
    output = []
    # replace protein ID with ID+gene name
    for line in lines:
        new_line = line
        for key in gene_names_dict.keys():
            new_line = new_line.replace(key, gene_names_dict[key])
            output.append(new_line)

# gos=[]
# for line in lines:
# 	for key in GOs_dict.keys():
# 		if key in line:
# 			gos.append(GOs_dict[key])

with open(gene_trees_appended, "w") as new_trees:
    for line in output:
        new_trees.write(line)
