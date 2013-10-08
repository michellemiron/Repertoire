#!/bin/bash

prefix=( TCR60 )

indir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/samples # directory with input files
outdir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/trimmed
mkdir -p $outdir  # directory to place split files

# split into pieces
for i in ${prefix[@]}
do
	echo "split -l 2000000 $indir/$i*R1* $outdir/tmp/$i.R1.part." | qsub -N $i.R1 -cwd -l mem=5G,time=3:: -o out -e err
	echo "split -l 2000000 $indir/$i*R2* $outdir/tmp/$i.R2.part." | qsub -N $i.R2 -cwd -l mem=5G,time=3:: -o out -e err
done
