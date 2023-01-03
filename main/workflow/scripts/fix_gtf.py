#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Natasha Picciani
Oct 31, 2022
usage: fix_gtf.py [transcriptome] [old_gff]

"""

import sys
import re
from Bio import SeqIO
import pandas as pd
import os

transcriptome_name = sys.argv[1]
gfffile = sys.argv[2]

outfile_tmp = f"{gfffile}.fixed.tmp.gff"
outfile = f"{gfffile}.fixed.gff"


def get_length(transcript):
    return str(transcript_length_dict[transcript])


# Create dictionary with transcript id:complete transcript name
# PB.1.1(key): PB.1.1|100013_0|path0:0-1954(+)|transcript/123808 (value)
transcripts = {}
transcript_length_dict = {}
with open(transcriptome_name) as transcriptome:
    for record in SeqIO.parse(transcriptome, "fasta"):
        transcript_shortid = record.id.split("|")[0]
        transcripts[transcript_shortid] = record.id
        transcript_length_dict[record.id] = len(record)

# Fix contig name - make it complete like in original transcriptome file
with open(outfile_tmp, "w") as gff_tmp:
    with open(gfffile) as gff:
        for line in gff:
            columns = line.split("\t")
            transcript_id = re.search('transcript_id "(PB.\d+.\d+)";', columns[8])
            transcript_complete = transcripts.get(transcript_id.group(1))
            columns[0] = transcript_complete
            if columns[2] == "exon":  # include lines where column 2 (third) is "exon"
                gff_tmp.write("\t".join(columns))

with open(outfile_tmp) as gff_tmp:
    gff_table = pd.read_table(gff_tmp, header=None)
    gff_table[9] = 1  # add column with value 1 in all rows
    gff_table[10] = gff_table[0].apply(get_length)
    newgff = gff_table[[0, 1, 2, 9, 10, 5, 6, 7, 8]]
    newgff.drop_duplicates(inplace=True)
    newgff.to_csv(outfile, sep="\t", quoting=3, header=None, index=False)

# os.remove(outfile_tmp)
