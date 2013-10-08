#!/bin/bash

kmerdir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/kmerfiles

:<< 'END'
for i in $kmerdir/*.kmers
do
	r=`echo $(basename $i .kmers).sorted` # | cut -f1 -d'_'`.sorted
	qsub -N `basename $i` -o log/`basename $i`.kcol -e log/`basename $i`.kcol splitsort_kmers.sh $kmerdir $i $r
	#echo "sort $i -T $kmerdir | uniq -c | sort -nrk1 -T $kmerdir > $i.frequencies" # | qsub -N `basename $i` -cwd -l mem=20G,time=48:: -o log/kmercol.out -e log/kmercol.err
done
END

#:<< 'END'
echo "sort -m $kmerdir/*.sorted/* -T $kmerdir > $kmerdir/humankmers.complete" | qsub -cwd -l mem=5G,time=10:: -o log/kmer.merge.out -e log/kmer.merge.err
# uniq -c > $i.counts" | qsub -N `basename $i` -cwd -l mem=5G,time=10:: -o log/kmer.unique.out -e log/kmer.unique.err
#END

