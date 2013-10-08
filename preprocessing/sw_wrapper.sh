prefix=( TCR22 TCR23 )
#prefix=( TCR22 TCR23 TCR24 TCR25 TCR26 TCR27 TCR28 TCR29 TCR30 TCR31 TCR32 TCR33 TCR34 TCR35 TCR36 TCR37 TCR38 TCR40 TCR41 TCR42 TCR43 ) #TCRX
#prefix=( TCR40 TCR41 TCR42 TCR43 )

#:<<'END'
# split and perform sw alignment
indir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/trimmed # directory with input files
#mkdir -p /ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/trimmed/tmp
for i in ${prefix[@]}
do
	fs=`ls $indir/tmp | grep $i.R1.part. | cut -f4 -d'.'`
	for j in $fs
	do
	   f1=$indir/tmp/$i.R1.part.$j
	   f2=$indir/tmp/$i.R2.part.$j
	  echo $f1
	   echo "./runsw 4 6 $f1 $f2 250 45 ${i}.${j}" | qsub -cwd -N $i.$j -l mem=5G,time=2:: -o /ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/preprocessing/log -e /ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/preprocessing/log
	done
done
#END

:<<'END'
# merge and cleanup
datadir="/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/errorcheck2"
for i in ${prefix[@]}
do
	cat $datadir/$i*.matches > $datadir/${i}_full.matches
	cat $datadir/$i*.mismatches > $datadir/${i}_full.mismatches
	cat $datadir/$i*.discard > $datadir/${i}_full.discard
	#rm $datadir/$i.*.matches
	#rm $datadir/$i.*.mismatches
	#rm $datadir/$i.*.discard
done
END
