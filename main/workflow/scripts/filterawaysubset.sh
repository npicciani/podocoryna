#!/bin/bash

filter_away_subset.py hq_transcripts.fasta.no5merge.collapsed
filter_away_subset.py Male_Med.Cogent.hq_transcripts.fasta.no5merge.collapsed
grep -c ">" Male_Med.Cogent.hq_transcripts.fasta.no5merge.collapsed.rep.fa
grep -c ">" Male_Med.Cogent.hq_transcripts.fasta.no5merge.collapsed.filtered.rep.fa
