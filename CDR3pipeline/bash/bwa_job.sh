#!/bin/bash
#$ -cwd
#$ -o log
#$ -e log


if [ "$#" != "5" ]; then
        echo 'usage: ' $0 '  <data file>  <reference>  <output directory> <BWA executable> <SAMTOOLS executable>'
        exit
fi

time1=$( date "+%s" )

echo [BEGIN] `date`
echo [MACHINE] `hostname`

#BWA=/ifs/scratch/c2b2/ys_lab/aps2157/software/bwa
#SAMTOOLS=/ifs/scratch/c2b2/ys_lab/bg2178/shared/software/samtools-1.0/samtools

fqfile=$1
ref=$2
outdir=$3
BWA=$4
SAMTOOLS=$5

name=${fqfile##*/}

if [ ! -d "$outdir" ]; then
	mkdir $outdir
fi

$BWA mem -t 4 -r 0.1 -E 1 -O 4 -B 3 -w 100 $ref $fqfile > $outdir/$name.sam
$SAMTOOLS view -bS $outdir/$name.sam > $outdir/$name.bam
rm $outdir/$name.sam

time2=$( date "+%s" )
echo [deltat] $(( $time2 - $time1 ))

echo [END] `date`
