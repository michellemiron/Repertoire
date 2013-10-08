vjdir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/human/VJcounts # $1
seg=$1 # V or J
type=$2 # all or 1000

# human B chain

# blood, blood+brain, GBM, LGG
blood=( TCR9 TCR27 TCR15 TCR13 TCR7 TCR23 TCR25 TCR29 )
sample_type_bl=( blood blood+brain GBM1 GBM2 LGG1 LGG2 LGG3 LGG4 )
brain=( TCR26 TCR14 TCR12 TCR6 TCR22 TCR24 TCR28 )
sample_type_br=( blood+brain GBM1 GBM2 LGG1 LGG2 LGG3 LGG4  )

f=`echo $vjdir/${blood[$1]}*Bchain*.$2.${seg}use.txt`
header=`cut -f1 $f | tr "\n" "\t"`

echo -e "name\ttissue\tsample\t$header">humanBchain.$2.${seg}.tsv

n=$[${#blood[@]}-1]
for i in $(seq 0 $n)
do
	f=`echo $vjdir/${blood[$i]}*Bchain*.$2.${seg}use.txt`
	vals=`cat $f | cut -f2 | tr "\n" "\t"`
	echo -e "${blood[$i]}\tblood\t${sample_type_bl[$i]}\t$vals" >> humanBchain.$2.${seg}.tsv
done 

n=$[${#brain[@]}-1]
for i in $(seq 0 $n)
do
        f=`echo $vjdir/${brain[$i]}*Bchain*.$2.${seg}use.txt`
        vals=`cat $f | cut -f2 | tr "\n" "\t"`
        echo -e "${brain[$i]}\tbrain\t${sample_type_br[$i]}\t$vals" >> humanBchain.$2.${seg}.tsv
done


# human A chain

# blood, blood+brain, GBM, LGG
blood=( TCR8 TCR26 TCR14 TCR12 TCR6 TCR22 TCR24 TCR28 ) 
sample_type_bl=( blood blood+brain GBM1 GBM2 LGG1 LGG2 LGG3 LGG4 )
brain=( TCR27 TCR15 TCR13 TCR7 TCR23 TCR25 TCR29 )
sample_type_br=( blood+brain GBM1 GBM2 LGG1 LGG2 LGG3 LGG4  )

f=`echo $vjdir/${blood[$1]}*Achain*.$2.${seg}use.txt`
header=`cut -f1 $f | tr "\n" "\t"`


echo -e "name\ttissue\tsample\t$header">humanAchain.$2.${seg}.tsv

n=$[${#blood[@]}-1]
for i in $(seq 0 $n)
do
        f=`echo $vjdir/${blood[$i]}*Achain*.$2.${seg}use.txt`
        vals=`cat $f | cut -f2 | tr "\n" "\t"`
        echo -e "${blood[$i]}\tblood\t${sample_type_bl[$i]}\t$vals" >> humanAchain.$2.${seg}.tsv
done

n=$[${#brain[@]}-1]
for i in $(seq 0 $n)
do
        f=`echo $vjdir/${brain[$i]}*Achain*.$2.${seg}use.txt`
        vals=`cat $f | cut -f2 | tr "\n" "\t"`
        echo -e "${brain[$i]}\tbrain\t${sample_type_br[$i]}\t$vals" >> humanAchain.$2.${seg}.tsv
done
