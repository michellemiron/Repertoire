#!/bin/bash
#$ -cwd
#$ -l mem=10G,time=48::
#$ -o log
#$ -e log

# Main pipeline to run the alpha chain analysis


if [ "$#" != "4" ]; then
        echo 'usage: ' $0 '  <data file> <output directory> <SAMTOOLS executable> <readfile>'
        exit
fi


time1=$( date "+%s" )

echo [BEGIN] `date`
echo [MACHINE] `hostname`

sample=$1 # sample name (bam file)
TRAoutdir=$2/`basename $1`.Achain # output directory


#SAMTOOLS=/ifs/scratch/c2b2/ys_lab/bg2178/shared/software/samtools-1.0/samtools
SAMTOOLS=$3
mkdir -p $TRAoutdir

$SAMTOOLS view -F 4 $sample > $TRAoutdir/`basename $sample .bam`.sam
python ./code/doTR/modify_sam.py $TRAoutdir/`basename $sample .bam`.sam $4 $TRAoutdir
perl ./code/doTR/ReadSam.CDR3.TRA.pl $TRAoutdir `basename $sample .bam`.sam  2>temp  #find VJ cassette and CDR3 region
perl ./code/doTR/ReadCDR3.VJ.TRA.pl $TRAoutdir `basename $sample .bam`.sam.CDR3.VJ.seq 2>temp  #Get reading frame and translate
perl ./code/doTR/CategorizeCDR3.TRA.pl $TRAoutdir `basename $sample .bam`.sam.CDR3.VJ.seq.prot.txt  2>temp # collect output
#rm $TRAoutdir/`basename $sample .bam`.sam

time2=$( date "+%s" )
echo [deltat] $(( $time2 - $time1 ))

echo [END] `date`
