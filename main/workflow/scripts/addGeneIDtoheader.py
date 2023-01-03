#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Natasha Picciani
last modified on May 02, 2022

Add partition information from Cogent run to fasta file header
Usage: addGeneIDtoheader.py filename.fasta final.partition.txt outputDirectory
Arguments:

	filename.fasta -- fasta file with isoform nucleotide sequences
	final.partition.txt -- output partition file from Cogent run with partition, size and isoform members
	outputDirectory -- path to output directory

"""
import pandas as pd
import re
import sys
from pathlib import PurePosixPath
from Bio import SeqIO

# User inputs
isoSeqFile = sys.argv[1]  # isoseq file with isoforms
partitionFile = sys.argv[2]  # Cogent partition file
outDir = sys.argv[3]  # output directory

## Read partition file and convert to a dictionary with isoforms as keys, partitions as values
partitions = pd.read_table(
    partitionFile, header=0
)  # read partition file as a dataframe without headers
partitions = partitions.drop(columns="Size")
partitionsDic = partitions.set_index("Partition").to_dict()["Members"]

q = dict()
isoDic = dict()

for key, value in partitionsDic.items():
    l = str(value).split(",")
    q[key] = l
    for iso in l:
        isoDic[iso] = key

## Export file with partition assignments per isoform
outfile = outDir + "/" + "final.partition.per.isoform.txt"
with open(outfile, "w") as out:
    for key, value in isoDic.items():
        out.write(key + "\t" + value + "\n")

## Export new fastafile with partition numbers in the headers. Repeat isoform name as partition assignment for isoforms without partition
outfile2 = outDir + "/" + PurePosixPath(isoSeqFile).name + ".fixed"
with open(isoSeqFile) as original, open(outfile2, "w") as corrected:
    for record in SeqIO.parse(original, "fasta"):
        if record.id in isoDic.keys():
            geneID = isoDic[record.id]
            record.id = record.id + "-Partition_" + geneID
        else:
            record.id = record.id + "-Partition_" + record.id
        SeqIO.write(record, corrected, "fasta")
