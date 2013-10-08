vjdir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/human/VJcounts # $1
seg=$1 # V or J
type=$2 # all or 1000

# human B chain

# blood, blood+brain, GBM, LGG
blood=( TCR40 TCR41 TCR42 TCR43 TCR7 TCR6 )
sample_type_bl=( NKT CD8 Treg CD4 LGG-blood LGG-brain)
brain=( SIMS0101-CD4-CD8 SIMS0101-CD8+CD56+ SIMS0101-CD4+CD25+ SIMS0101-PBMC SIMS0101-TUMOR )
sample_type_br=( NKT CD8 Treg LGG-blood LGG-brain )

f=`echo $vjdir/${blood[$1]}*Bchain*.$2.${seg}use.txt`
header=`cut -f1 $f | tr "\n" "\t"`

echo -e "name\ttissue\tsample\t$header">TCR40-43.$2.${seg}.tsv

n=$[${#blood[@]}-1]
for i in $(seq 0 $n)
do
	f=`echo $vjdir/${blood[$i]}*Bchain*.$2.${seg}use.txt`
	vals=`cat $f | cut -f2 | tr "\n" "\t"`
	echo -e "${blood[$i]}\tblood\t${sample_type_bl[$i]}\t$vals" >> TCR40-43.$2.${seg}.tsv
done 

n=$[${#brain[@]}-1]
for i in $(seq 0 $n)
do
        f=`echo $vjdir/${brain[$i]}*.$2.${seg}use.txt`
        vals=`cat $f | cut -f2 | tr "\n" "\t"`
        echo -e "${brain[$i]}\tbrain\t${sample_type_br[$i]}\t$vals" >> TCR40-43.$2.${seg}.tsv
done
