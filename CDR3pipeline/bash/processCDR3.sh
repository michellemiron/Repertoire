#!/bin/bash

if [ "$#" != "5" ]; then
        echo 'usage: ' $0 '  <data name>  <data directory>  <output directory>  <chain: A for alpha, B for beta> <awk>'
        exit
fi

prefix=$1
datadir=$2
outdir=$3
chain=$4
aw=$5
mkdir -p $outdir

./bash/processCDR3s_job.sh $prefix $datadir $outdir $chain $aw
