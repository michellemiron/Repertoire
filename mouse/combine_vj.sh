vjdir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/mouse/VJcounts # $1
seg=$1 # V or J
type=$2 # all or 1000

# human B chain

# blood, blood+brain, Day 21, Day 28, Day 42, End stage
blood=( TCR11 TCR17 TCR31 TCR19 TCR21 TCR33 TCR39 )
sample_type_bl=( blood Day21-1 Day21-2 Day28 Day42-1 Day42-2 Endstage  )
brain=( TCR16 TCR30 TCR18 TCR20 TCR32 TCR38 )
sample_type_br=( Day21-1 Day21-2 Day28 Day42-1 Day42-2 Endstage )

f=`echo $vjdir/${blood[$1]}*Bchain*.$2.${seg}use.txt`
header=`cut -f1 $f | tr "\n" "\t"`

echo -e "name\ttissue\tsample\t$header">mouseBchain.$2.${seg}.tsv

n=$[${#blood[@]}-1]
for i in $(seq 0 $n)
do
	f=`echo $vjdir/${blood[$i]}*Bchain*.$2.${seg}use.txt`
	vals=`cat $f | cut -f2 | tr "\n" "\t"`
	echo -e "${blood[$i]}\tblood\t${sample_type_bl[$i]}\t$vals" >> mouseBchain.$2.${seg}.tsv
done 

n=$[${#brain[@]}-1]
for i in $(seq 0 $n)
do
        f=`echo $vjdir/${brain[$i]}*Bchain*.$2.${seg}use.txt`
        vals=`cat $f | cut -f2 | tr "\n" "\t"`
        echo -e "${brain[$i]}\tbrain\t${sample_type_br[$i]}\t$vals" >> mouseBchain.$2.${seg}.tsv
done


# human A chain

# blood, blood+brain, Day 21, Day 28, Day 42, End stage
blood=( TCR10 TCR16 TCR30 TCR18 TCR20 TCR32 TCR38 )
sample_type_bl=( blood Day21-1 Day21-2 Day28 Day42-1 Day42-2 Endstage  )
brain=( TCR17 TCR31 TCR19 TCR21 TCR33 TCR39 )
sample_type_br=( Day21-1 Day21-2 Day28 Day42-1 Day42-2 Endstage )

f=`echo $vjdir/${blood[$1]}*Achain*.$2.${seg}use.txt`
header=`cut -f1 $f | tr "\n" "\t"`


echo -e "name\ttissue\tsample\t$header">mouseAchain.$2.${seg}.tsv

n=$[${#blood[@]}-1]
for i in $(seq 0 $n)
do
        f=`echo $vjdir/${blood[$i]}*Achain*.$2.${seg}use.txt`
        vals=`cat $f | cut -f2 | tr "\n" "\t"`
        echo -e "${blood[$i]}\tblood\t${sample_type_bl[$i]}\t$vals" >> mouseAchain.$2.${seg}.tsv
done

n=$[${#brain[@]}-1]
for i in $(seq 0 $n)
do
        f=`echo $vjdir/${brain[$i]}*Achain*.$2.${seg}use.txt`
        vals=`cat $f | cut -f2 | tr "\n" "\t"`
        echo -e "${brain[$i]}\tbrain\t${sample_type_br[$i]}\t$vals" >> mouseAchain.$2.${seg}.tsv
done
