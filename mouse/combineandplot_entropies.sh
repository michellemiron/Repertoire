entropydir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/mouse/entropy

# mouse B chain

# blood, blood+brain, Day 21, Day 28, Day 42, End stage
blood=( TCR11 TCR17 TCR31 TCR19 TCR21 TCR33 TCR39 )
sample_type_bl=( blood Day21-1 Day21-2 Day28 Day42-1 Day42-2 Endstage  )
brain=( TCR16 TCR30 TCR18 TCR20 TCR32 TCR38 )
sample_type_br=( Day21-1 Day21-2 Day28 Day42-1 Day42-2 Endstage )

echo -e "name\ttissue\tsample\ttotal_reads\tentropy_CDR3_total\tentropy_CDR3_1000\tentropy_VJ_total\tentropy_VJ_1000\tentropy_V_total\tentropy_V_1000\tentropy_J_total\tentropy_J_1000">mouseBchain.entropies.tsv

n=$[${#blood[@]}-1]
for i in $(seq 0 $n)
do
	reads=`cat  $entropydir/${blood[$i]}*.Bchain*.all.Hcdr3 | cut -f1`
	Hcdr3=`cat  $entropydir/${blood[$i]}*.Bchain*.all.Hcdr3 | cut -f2`
	Hcdr3_1000=`cat $entropydir/${blood[$i]}*.Bchain*.1000.Hcdr3 | cut -f2`
	Hvj=`cat $entropydir/${blood[$i]}*.Bchain*.all.Hvj | cut -f2`
	Hvj_1000=`cat $entropydir/${blood[$i]}*.Bchain*.1000.Hvj | cut -f2`  
        Hv=`cat $entropydir/${blood[$i]}*.Bchain*.all.Hv | cut -f2`          
        Hv_1000=`cat $entropydir/${blood[$i]}*.Bchain*.1000.Hv | cut -f2`
        Hj=`cat $entropydir/${blood[$i]}*.Bchain*.all.Hj | cut -f2`          
        Hj_1000=`cat $entropydir/${blood[$i]}*.Bchain*.1000.Hj | cut -f2`
	echo -e "${blood[$i]}\tblood\t${sample_type_bl[$i]}\t$reads\t$Hcdr3\t$Hcdr3_1000\t$Hvj\t$Hvj_1000\t$Hv\t$Hv_1000\t$Hj\t$Hj_1000" >> mouseBchain.entropies.tsv
done

n=$[${#brain[@]}-1]
for i in $(seq 0 $n)
do
	reads=`cat $entropydir/${brain[$i]}*.Bchain*.all.Hcdr3 | cut -f1`
        Hcdr3=`cat $entropydir/${brain[$i]}*.Bchain*.all.Hcdr3 | cut -f2`
        Hcdr3_1000=`cat $entropydir/${brain[$i]}*.Bchain*.1000.Hcdr3 | cut -f2`
        Hvj=`cat $entropydir/${brain[$i]}*.Bchain*.all.Hvj | cut -f2`
        Hvj_1000=`cat $entropydir/${brain[$i]}*.Bchain*.1000.Hvj | cut -f2`
        Hv=`cat $entropydir/${brain[$i]}*.Bchain*.all.Hv | cut -f2`
        Hv_1000=`cat $entropydir/${brain[$i]}*.Bchain*.1000.Hv | cut -f2`
        Hj=`cat $entropydir/${brain[$i]}*.Bchain*.all.Hj | cut -f2`
        Hj_1000=`cat $entropydir/${brain[$i]}*.Bchain*.1000.Hj | cut -f2`
	echo -e "${brain[$i]}\tbrain\t${sample_type_br[$i]}\t$reads\t$Hcdr3\t$Hcdr3_1000\t$Hvj\t$Hvj_1000\t$Hv\t$Hv_1000\t$Hj\t$Hj_1000" >> mouseBchain.entropies.tsv
done

# human A chain

# blood, blood+brain, Day 21, Day 28, Day 42, End stage
blood=( TCR10 TCR16 TCR30 TCR18 TCR20 TCR32 TCR38 )
sample_type_bl=( blood Day21-1 Day21-2 Day28 Day42-1 Day42-2 Endstage  )
brain=( TCR17 TCR31 TCR19 TCR21 TCR33 TCR39 )
sample_type_br=( Day21-1 Day21-2 Day28 Day42-1 Day42-2 Endstage )

echo -e "name\ttissue\tsample\ttotal_reads\tentropy_CDR3_total\tentropy_CDR3_1000\tentropy_VJ_total\tentropy_VJ_1000\tentropy_V_total\tentropy_V_1000\tentropy_J_total\tentropy_J_1000">mouseAchain.entropies.tsv

n=$[${#blood[@]}-1]
for i in $(seq 0 $n)
do
        reads=`cat  $entropydir/${blood[$i]}*.Achain*.all.Hcdr3 | cut -f1`
        Hcdr3=`cat  $entropydir/${blood[$i]}*.Achain*.all.Hcdr3 | cut -f2`
        Hcdr3_1000=`cat $entropydir/${blood[$i]}*.Achain*.1000.Hcdr3 | cut -f2`
        Hvj=`cat $entropydir/${blood[$i]}*.Achain*.all.Hvj | cut -f2`
        Hvj_1000=`cat $entropydir/${blood[$i]}*.Achain*.1000.Hvj | cut -f2`
        Hv=`cat $entropydir/${blood[$i]}*.Achain*.all.Hv | cut -f2`
        Hv_1000=`cat $entropydir/${blood[$i]}*.Achain*.1000.Hv | cut -f2`
        Hj=`cat $entropydir/${blood[$i]}*.Achain*.all.Hj | cut -f2`
        Hj_1000=`cat $entropydir/${blood[$i]}*.Achain*.1000.Hj | cut -f2`
        echo -e "${blood[$i]}\tblood\t${sample_type_bl[$i]}\t$reads\t$Hcdr3\t$Hcdr3_1000\t$Hvj\t$Hvj_1000\t$Hv\t$Hv_1000\t$Hj\t$Hj_1000" >> mouseAchain.entropies.tsv
done

n=$[${#brain[@]}-1]
for i in $(seq 0 $n)
do
        reads=`cat $entropydir/${brain[$i]}*.Achain*.all.Hcdr3 | cut -f1`
        Hcdr3=`cat $entropydir/${brain[$i]}*.Achain*.all.Hcdr3 | cut -f2`
        Hcdr3_1000=`cat $entropydir/${brain[$i]}*.Achain*.1000.Hcdr3 | cut -f2`
        Hvj=`cat $entropydir/${brain[$i]}*.Achain*.all.Hvj | cut -f2`
        Hvj_1000=`cat $entropydir/${brain[$i]}*.Achain*.1000.Hvj | cut -f2`
        Hv=`cat $entropydir/${brain[$i]}*.Achain*.all.Hv | cut -f2`
        Hv_1000=`cat $entropydir/${brain[$i]}*.Achain*.1000.Hv | cut -f2`
        Hj=`cat $entropydir/${brain[$i]}*.Achain*.all.Hj | cut -f2`
        Hj_1000=`cat $entropydir/${brain[$i]}*.Achain*.1000.Hj | cut -f2`
        echo -e "${brain[$i]}\tbrain\t${sample_type_br[$i]}\t$reads\t$Hcdr3\t$Hcdr3_1000\t$Hvj\t$Hvj_1000\t$Hv\t$Hv_1000\t$Hj\t$Hj_1000" >> mouseAchain.entropies.tsv
done

