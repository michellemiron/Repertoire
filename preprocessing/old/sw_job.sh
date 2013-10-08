#!/bin/bash
#$ -cwd
#$ -o log
#$ -e log

time1=$( date "+%s" )

echo [BEGIN] `date`
echo [MACHINE] `hostname`

splitfile=$1
outdir=$2

outfile=$outdir/$(basename $splitfile| cut -f1 -d'.')

echo $outfile
awk 'NR%4==2 { print $0 } ' OFS=\\t $splitfile | \
# awk 'NR%4==1 { header=$1 } NR%4==2 { print header , $0 } ' OFS=\\t $splitfile  | \
while read s1 s2
do
	./runsw 4 6 $s1 $(echo $s2 | tr ACGT TGCA | rev) 250 | awk '{print $1}' > $outfile
done
#else ./runsw 4 6 \$s1 \$(echo \$s2 | tr [ACGT] [TGCA] | rev) 250 | awk -v head=\"\$h\" -v outdir=\"$outdir\" -v id=\"$i\" '{if(\$1==\$2) print head\"\n\"\$1 >> outdir\"/\"id\".match\"; else print head\"\n\"\$0 >> outdir\"/\"id\".mismatch\"}'
#fi; rm $j done

time2=$( date "+%s" )
echo [deltat] $(( $time2 - $time1 ))

echo [END] `date`

