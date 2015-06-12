#!/bin/bash
#$ -cwd
#$ -l mem=10G,time=48::
#$ -o log
#$ -e log

if [ "$#" != "4" ]; then
        echo 'usage: ' $0 '  <data file>  <output directory> <SAMTOOLS executable> <readfile>'
        exit
fi


time1=$( date "+%s" )

echo [BEGIN] `date`
echo [MACHINE] `hostname`
                              
sample=$1 # sample name (bam file)
TRBoutdir=$2/`basename $1`.Bchain # output directory

#SAMTOOLS=/ifs/scratch/c2b2/ys_lab/bg2178/shared/software/samtools-1.0/samtools
SAMTOOLS=$3
mkdir -p $TRBoutdir

$SAMTOOLS view -F 4 $sample > $TRBoutdir/`basename $sample .bam`.sam
python ./code/doTR/modify_sam.py $TRBoutdir/`basename $sample .bam`.sam $4 $TRBoutdir
perl ./code/doTR/ReadSam.CDR3.TRB.pl $TRBoutdir `basename $sample .bam`.sam 2>temp
perl ./code/doTR/ReadCDR3.VJ.TRB.pl $TRBoutdir `basename $sample .bam`.sam.CDR3.VJ.seq 2>temp
perl ./code/doTR/CategorizeCDR3.TRB.pl $TRBoutdir `basename $sample .bam`.sam.CDR3.VJ.seq.prot.txt 2>temp
rm $TRBoutdir/`basename $sample .bam`.sam


time2=$( date "+%s" )
echo [deltat] $(( $time2 - $time1 ))

echo [END] `date`
