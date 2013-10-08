#!/bin/bash

prefix=$1
datadir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/bamfiles2
fastapath=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/fastas/human_g1k_v37.MaskTRB.TRBp8.fasta
outpath=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/fastqcorrected

fs=`ls $datadir | grep $prefix | cut -f4 -d"."| sort -u`
for i in $fs
do
	#for j in $datadir/$prefix*$i*
	#do
		f1=$datadir/$prefix.R1.part.$i.bam
		f2=$datadir/$prefix.R2.part.$i.bam
		echo $f1
		outname=$prefix.$i.corrected
		#sh errorcorrectionjob.sh $f1 $f2 $datadir $fastapath $outpath $outname
		qsub -N $prefix.$i errorcorrectionjob.sh $f1 $f2 $datadir $fastapath $outpath $outname
	#done
done

