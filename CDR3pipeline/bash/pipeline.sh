#!/bin/bash

if [ "$#" != "6" ]; then
        echo 'usage: ' $0 '  <data name>  <data directory>  <output directory>  <chain: A for alpha, B for beta> <SAMTOOLS executable> <Readfile>'
        exit
fi


prefix=$1
datadir=$2
outdir=$3
chain=$4
SAMTOOLS=$5
Readfile=$6
if [ $chain == A ]; then
	code=./bash/doTRA.sh
elif [ $chain == B ]; then
	code=./bash/doTRB.sh
else
	echo "Improper input chain = $chain"
fi


i=$datadir/$prefix


$code $i $outdir $SAMTOOLS $Readfile
