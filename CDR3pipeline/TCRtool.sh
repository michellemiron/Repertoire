#!/bin/bash

if [ $# -lt 2 ]; then
        echo "usage: ./TCRtool <sample file> <reference genome> -options"
	echo "options:"
	echo "         -o   Output directory."
	echo "              [default]./TCRoutput"
	echo "         -C   The DNA chain needed to be analysed. Can be 'A', 'B' or 'both'."
	echo "              [default]both"
	echo "         -B   BWA executable's path."
	echo "              [default]BWA"
	echo "         -S   samtools executable's path."
	echo "              [default]samtools"
	echo "         -A   awk GNU version executable's path."
	echo "              [default]awk"
        exit
fi


data=$1
chain='both'
odir='./TCRoutput'
ref=$2
bwa='BWA'
samtools='samtools'
aw='awk'


OPTIND=3
while getopts "o:B:S:A:C:" arg
do
        case $arg in
                o) odir=$OPTARG;;
                B) bwa=$OPTARG;;
                S) samtools=$OPTARG;;
		A) aw=$OPTARG;;
		C) chain=$OPTARG;;
                \?)        echo "usage: ./TCRtool <sample file> <reference genome> -options"
        echo "options:"
        echo "         -o   Output directory."
        echo "              [default]./TCRoutput"
        echo "         -C   The DNA chain needed to be analysed. Can be 'A', 'B' or 'both'."
        echo "              [default]both"
        echo "         -B   BWA executable's path."
        echo "              [default]BWA"
        echo "         -S   samtools executable's path."
        echo "              [default]samtools"
        echo "         -A   awk GNU version executable's path."
        echo "              [default]awk"

                    exit;;
        esac
done


echo 'Data: '$data
echo 'References Genome: '$ref
echo 'Output Directory: '$odir
echo 'Chain: '$chain
echo 'BWA executable: '$bwa
echo 'samtools executable: '$samtools
echo 'awk executable: ' $aw


mkdir -p $odir/createddirforTCRtools1
mkdir -p $odir/createddirforTCRtools2

if [ "$chain" == "both" ];then
chains='A'
echo '---------------------------------------------------------'
echo 'Reading Files'
echo '---------------------------------------------------------'
./bash/bwa_job.sh $data $ref $odir/createddirforTCRtools1 $bwa $samtools
echo '---------------------------------------------------------'
echo 'Processing alpha chain'
echo '---------------------------------------------------------'
./bash/pipeline.sh ${data##*/}.bam $odir/createddirforTCRtools1 $odir/createddirforTCRtools2 $chains $samtools $data
./bash/processCDR3.sh ${data##*/}.bam $odir/createddirforTCRtools2 $odir $chains $aw

chains='B'
echo '---------------------------------------------------------'
echo 'Processing beta chain'
echo '---------------------------------------------------------'
./bash/pipeline.sh ${data##*/}.bam $odir/createddirforTCRtools1 $odir/createddirforTCRtools2 $chains $samtools $data
./bash/processCDR3.sh ${data##*/}.bam $odir/createddirforTCRtools2 $odir $chains $aw

rm -rf $odir/createddirforTCRtools1
rm -rf $odir/createddirforTCRtools2
rm -f temp

echo '---------------------------------------------------------'
echo 'Finished'
echo '---------------------------------------------------------'
fi


if [ "$chain" == "A" ];then
chains='A'
echo '---------------------------------------------------------'
echo 'Reading Files'
echo '---------------------------------------------------------'
./bash/bwa_job.sh $data $ref $odir/createddirforTCRtools1 $bwa $samtools
echo '---------------------------------------------------------'
echo 'Processing alpha chain'
echo '---------------------------------------------------------'
./bash/pipeline.sh ${data##*/}.bam $odir/createddirforTCRtools1 $odir/createddirforTCRtools2 $chains $samtools $data
./bash/processCDR3.sh ${data##*/}.bam $odir/createddirforTCRtools2 $odir $chains $aw


rm -rf $odir/createddirforTCRtools1
rm -rf $odir/createddirforTCRtools2
rm -f temp

echo '---------------------------------------------------------'
echo 'Finished'
echo '---------------------------------------------------------'
fi


if [ "$chain" == "B" ];then
echo '---------------------------------------------------------'
echo 'Reading Files'
echo '---------------------------------------------------------'
./bash/bwa_job.sh $data $ref $odir/createddirforTCRtools1 $bwa $samtools
chains='B'
echo '---------------------------------------------------------'
echo 'Processing beta chain'
echo '---------------------------------------------------------'
./bash/pipeline.sh ${data##*/}.bam $odir/createddirforTCRtools1 $odir/createddirforTCRtools2 $chains $samtools $data
./bash/processCDR3.sh ${data##*/}.bam $odir/createddirforTCRtools2 $odir $chains $aw

rm -rf $odir/createddirforTCRtools1
rm -rf $odir/createddirforTCRtools2
rm -f temp

echo '---------------------------------------------------------'
echo 'Finished'
echo '---------------------------------------------------------'
fi

