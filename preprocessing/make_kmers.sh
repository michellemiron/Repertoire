#!/bin/bash

outdir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/kmerfiles

for i in `cat /ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/human/paired_end_samples`
do
	for j in /ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/trimmed/${i}_*.fastq
	do
		> $outdir/`basename $j .trimed.fastq`.kmers
		echo "cat $j | grep -v removed | awk '{if(NR%4==2) print \$0}' | python read2kmer.py $outdir/`basename $j .trimed.fastq`.kmers" | qsub -N `basename $i` -cwd -l mem=10G,time=48:: -e log/kmer.err -o log/kmer.out
	done
done

#echo "cat /ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/trimmed/TCR4_S1_L001_R1_001.trimed.fastq | grep -v removed| awk '{if(NR%4==2) print \$0}'  | python read2kmer.py test"  | qsub -N test -cwd -l mem=10G,time=24::

