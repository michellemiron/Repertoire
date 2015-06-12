#!/bin/bash
#$ -cwd
#$ -l mem=5G,time=2::
#$ -o log
#$ -e log


if [ "$#" != "5" ]; then
        echo 'usage: ' $0 '  <data name>  <data directory>  <output directory>  <chain: A for alpha, B for beta>' '<awk exe>'
        exit
fi



prefix=$1
datadir=$2
outdir=$3
chain=$4
aw=$5

# get all CDR3 sequences
i=$datadir/$prefix.${chain}chain
	cat $i/*.seq.prot.txt | cut -f4,5,12,15 | sort -k3 | uniq -c | $aw '{print $2"."$3"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"length($4)}' > $outdir/$prefix.${chain}chain.full 


# get productive CDR3 sequences
infile=$outdir/$prefix.${chain}chain.full
echo -e "VJCombine\tCopy\tVGeneName\tJGeneName\tCDR3\tSequence\tLengthOfCDR3" >$outdir/`basename $infile .full`.productive.tsv ; cat $infile |grep -v Out | grep -v "*" | $aw '{if(NF==7) print $0}' | sort -k2 -nr >> $outdir/`basename $infile .full`.productive.tsv;


# final error correction step
infile=$outdir/`basename $infile .full`.productive.tsv
Jmotif=./files/processCDR3s/human${chain}.Jmotif

echo -e "VJCombine\tCopy\tVGeneName\tJGeneName\tCDR3\tSequence" > $outdir/`basename $infile .productive.tsv`.final.tsv
sed 1d $infile |  python ./code/processCDR3s/fixNTs.py $Jmotif  |  $aw -f ./code/processCDR3s/merge.awk | tr @ \\t | sort -rnk6 | $aw '{print $1"\t"$6"\t"$2"\t"$3"\t"$4"\t"$5}' >> $outdir/`basename $infile .productive.tsv`.final.tsv

rm -f $outdir/$prefix.${chain}chain.full
rm -f $outdir/$prefix.${chain}chain.productive.tsv
mv $outdir/$prefix.${chain}chain.final.tsv $outdir/${prefix%.*}.${chain}chain.tsv

