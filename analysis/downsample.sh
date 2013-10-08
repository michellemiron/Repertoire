#!/bin/bash
#$ -cwd
#$ -o log
#$ -e log

source ~/.bashrc

time1=$( date "+%s" )

echo [BEGIN] `date`
echo [MACHINE] `hostname`

picardpath=/ifs/data/c2b2/ngs_lab/ngs/usr/picard-tools-1.65
datadir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/bamfiles

input=$1
fraction=$2
no=$3
seed=$[$time1+$RANDOM] # random int added to adjust seed for jobs started at roughly the same time

outdir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/bamfiles/downsampled/ds${fraction}
mkdir -p $outdir


java -Xmx5g -jar ${picardpath}/DownsampleSam.jar INPUT=$datadir/$input OUTPUT=${outdir}/`basename $input .bam`_${fraction}_$no.bam PROBABILITY=$fraction VALIDATION_STRINGENCY=SILENT RANDOM_SEED=$seed

time2=$( date "+%s" )
echo [deltat] $(( $time2 - $time1 ))

echo [END] `date`
