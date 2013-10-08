#!/bin/sh

processedfile=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/fastqcorrected/processedbams/TCRB/TCR62full.VJ.seq.prot.txt
processedbarcodes=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/fastqcorrected/processedbams/TCRB/TCR62.barcodes.txt
barcodelist=barcodes.tsv
motif=CCCA
outdir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/fastqcorrected/processedbams/TCRB

python demultiplex.py $processedfile $processedbarcodes $barcodelist $motif $outdir
