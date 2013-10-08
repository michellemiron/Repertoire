#!/bin/bash
#$ -cwd
#$ -l mem=10G,time=48::


time1=$( date "+%s" )

echo [BEGIN] `date`
echo [MACHINE] `hostname`

kmerdir=$1
fpath=$2
s=$3

mkdir -p $kmerdir/$s

split -l 249660952 $fpath $kmerdir/$s/`basename $s .sorted`

for i in $kmerdir/$s/*
do
	echo "sort -T $kmerdir/$s $i> $i.sorted ; rm $i" | qsub -N `basename $i` -cwd -l mem=8G,time=10:: -o log/kmersort.out -e log/kmersort.err
done

time2=$( date "+%s" )
echo [deltat] $(( $time2 - $time1 ))

echo [END] `date`
