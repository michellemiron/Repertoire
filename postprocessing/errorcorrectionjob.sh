#!bin/bash
#$ -l mem=5G,time=24:: 
#$ -cwd
#$ -o log
#$ -e log

time1=$( date "+%s" )
echo [BEGIN] `date`
echo [MACHINE] `hostname`
        
	
f1=$1
f2=$2
datadir=$3
fastapath=$4
outpath=$5
outname=$6
	samtools view $f1 > `basename $f1 .bam`.sam
        samtools view $f2 > `basename $f2 .bam`.sam
        ./mergereads `basename $f1 .bam`.sam /ifs/home/c2b2/ys_lab/bg2178/tcr-src/postprocessing 
        ./mergereads `basename $f2 .bam`.sam /ifs/home/c2b2/ys_lab/bg2178/tcr-src/postprocessing 
#echo $prefix.R1.$i.sam.merged  
        ./fixreads `basename $f1 .bam`.sam.merged `basename $f2 .bam`.sam.merged $fastapath $outpath $outname
	rm `basename $f1 .bam`.sam
	rm `basename $f2 .bam`.sam
	rm `basename $f1 .bam`.sam.merged
	rm `basename $f2 .bam`.sam.merged

time2=$( date "+%s" )
echo [deltat] $(( $time2 - $time1 ))

echo [END] `date`

