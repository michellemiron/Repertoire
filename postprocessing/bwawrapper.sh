#!/bin/bash

fastadir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/fastas
humanref=human_g1k_v37.MaskTRB.TRBp8.fasta
mouseref=mm10.MaskTraTrb.TRAp1TRBp1.fa

#mapreadydir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/errorcheck2
#mapreadydir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/immunoseq
#outdir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/errorcheck2/matchbams2
#outdir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/immunoseq/immunobams
mapreadydir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/fastqcorrected
outdir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/fastqcorrected/bamfiles
mkdir -p $outdir

# change these depending on samples
ref=$fastadir/$humanref
#fileid=( TCR22 TCR23 )
fileid=( TCR58 )

#ref=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/immunoseq/TCR7_10000seqs.fasta

for i in ${fileid[@]}
do
#fqfiles=$mapreadydir/TCR${i}.*.part*
fqfiles=$mapreadydir/$i*.fastq
#fqfiles=$mapreadydir/$i.be.*.fastq
for fq in $fqfiles
do
	echo $fq
	qsub -V -N bwa${i} -cwd -pe smp 4 -R y -l mem=5G,time=2:: bwa_job.sh $fq $ref $outdir
done
done
