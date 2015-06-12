#!/bin/bash

# Take a sample and split it up into piece of size 1000000 for quicker processing
if [ "$#" != "3" ]; then
	echo 'usage: ' $0 '  <sample name>  <sample directory>  <output directory>'
	exit
fi

sample=$1 # sample name
indir=$2 # sample directory
outdir=$3 # directory to place the output


mkdir -p $outdir  # make output directory if it doesn't exist

# split into pieces
echo "split -l 700000 $indir/$sample $outdir/`basename $sample .fastq`.part." | qsub -N $sample -cwd -o log -e log
