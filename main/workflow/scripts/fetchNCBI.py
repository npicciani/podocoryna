#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python 3.7
# -*- coding: utf-8 -*-

"""
created on June 3, 2019
by Natasha Picciani
usage: fetch.py list outputDir

This script pull sequences for a list of NCBI accession numbers.
Output is a file named "fastaRecords.fasta" and another with underscores
instead of spaces named "fastaRecords_underscores.fasta" in the output directory specified

"""
import sys
from Bio import Entrez
from urllib.error import HTTPError
import re


# Email settings for NCBI and user inputs
Entrez.email = "n_picciani@hotmail.com"  # user defined
ListPath = sys.argv[1]  # list of sequences to fetch
outDir = sys.argv[2]  # output directory

# Listing the IDs of sequences available on the NCBI nucleotide database from a list of taxa
accNum = []
with open(ListPath, "r") as seqList:
    for line in seqList:
        line = line.strip("\n")
        accNum.append(line)

# Pulling fasta files from list of accession numbers
fastaRecordsPath = outDir + "/fastaRecords.fasta"
fastaRecords = []
for item in accNum:
    # try/except block to prevent HTTPError /depengli
    try:
        handle = Entrez.efetch(
            db="protein", id=item, rettype="fasta", retmode="text"
        )  # automatically uses an HTTP POST if there are over 200 identifiers
    except HTTPError:
        time.sleep(5)
        handle = Entrez.efetch(db="protein", id=item, rettype="fasta", retmode="text")
    fastaRecords.extend(handle)
    handle.close()

with open(fastaRecordsPath, "w") as f:
    for item in fastaRecords:
        f.write("%s" % item)

print("Exported fasta file")

# Replacing white spaces in the fasta headers
pattern = re.compile(r" ")
replacement = r"_"
fastaFileUnderscore = outDir + "/fastaRecords_underscores.fasta"
outFile = open(fastaFileUnderscore, "w")

with open(fastaRecordsPath, "r") as file:
    for line in file:
        if line.startswith(">"):
            newline = pattern.sub(replacement, line)
            outFile.write(newline)
        else:
            outFile.write(line)

outFile.close()
