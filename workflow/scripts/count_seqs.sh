#!/bin/bash

FASTA=Podocoryna_carnea.collapsed.fasta
SCRIPTS=/home/nnp9/project/podocoryna/workflow/scripts/

cd /home/nnp9/scratch60/podocoryna/results/original_precogent/treeinform/
python $SCRIPTS/identify_transcripts.py ./threshold_0.5/$FASTA >> threshold_counts.txt
python $SCRIPTS/identify_transcripts.py ./threshold_1/$FASTA >> threshold_counts.txt
python $SCRIPTS/identify_transcripts.py ./threshold_2/$FASTA >> threshold_counts.txt
python $SCRIPTS/identify_transcripts.py ./threshold_3/$FASTA >> threshold_counts.txt
python $SCRIPTS/identify_transcripts.py ./threshold_4/$FASTA >> threshold_counts.txt
python $SCRIPTS/identify_transcripts.py ./threshold_5/$FASTA >> threshold_counts.txt
