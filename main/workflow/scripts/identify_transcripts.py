# -*- coding: utf-8 -*-
# Written by Natasha Picciani
# Oct 31, 2022
# usage: identify_transcripts.py fastafile

from Bio import SeqIO
import re
import sys

collapsed = sys.argv[1]
protein_seqs = []
transcript_seqs = set()
with open(collapsed) as fastafile:
    for record in SeqIO.parse(fastafile, "fasta"):
        protein_seqs.append(record.id)
        # transcript_id = re.search(r"(PB.\d+.\d+|)", record.id)        transcript_id = re.search(r"(transcript_\d+).p\d+ ", record.id)
        # transcript_seqs.add(transcript_id[0])
        transcript_id = re.search(r"(transcript_\d+).p\d+", record.id)
        transcript_seqs.add(transcript_id.group(1))
print(str(len(protein_seqs)) + "	" + str(len(transcript_seqs)))
