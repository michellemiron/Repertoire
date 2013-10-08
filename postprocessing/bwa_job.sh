#!/bin/bash
#$ -cwd
#$ -o log
#$ -e log

time1=$( date "+%s" )

echo [BEGIN] `date`
echo [MACHINE] `hostname`

BWA=/ifs/scratch/c2b2/ys_lab/aps2157/software/bwa

fqfile=$1
ref=$2
outdir=$3
echo $ref
echo $fqfile

name=`basename $fqfile .fastq`
$BWA bwasw -t 4 -z 3 -c 3 -r 1 -q 4 -w 100 $ref $fqfile > $outdir/$name.sam
samtools view -bS $outdir/$name.sam > $outdir/$name.bam
#rm $outdir/$name.sam


#TRA
#samtools view $a.bam |awk '{if ($3==14&&($4>=22090057||$4<=23021075)) print $0}' >$a.TRAhits.sam 

#TRB
#samtools view $a.bam |awk '{if ($3=="gi|114841177:91557-667340") print $0}' >$a.TRBhits.sam 


time2=$( date "+%s" )
echo [deltat] $(( $time2 - $time1 ))

echo [END] `date`
