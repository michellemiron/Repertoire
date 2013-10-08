entropydir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/human/entropy # $1

# iRep vs iSeq

# blood, blood+brain, GBM, LGG
blood=( TCR40 TCR41 TCR42 TCR43 TCR7 TCR6 )
sample_type_bl=( NKT CD8 Treg CD4 LGG-blood LGG-brain)
brain=( SIMS0101-CD4-CD8 SIMS0101-CD8+CD56+ SIMS0101-CD4+CD25+ SIMS0101-PBMC SIMS0101-TUMOR )
sample_type_br=( NKT CD8 Treg LGG-blood LGG-brain )

echo -e "name\ttissue\tsample\ttotal_reads\tentropy_CDR3_total\tentropy_CDR3_1000\tentropy_VJ_total\tentropy_VJ_1000\tentropy_V_total\tentropy_V_1000\tentropy_J_total\tentropy_J_1000">TCR40-43.entropies.tsv

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
	echo -e "${blood[$i]}\tblood\t${sample_type_bl[$i]}\t$reads\t$Hcdr3\t$Hcdr3_1000\t$Hvj\t$Hvj_1000\t$Hv\t$Hv_1000\t$Hj\t$Hj_1000" >> TCR40-43.entropies.tsv
done 

n=$[${#brain[@]}-1]
for i in $(seq 0 $n)
do
	reads=`cat $entropydir/${brain[$i]}*.all.Hcdr3 | cut -f1`
        Hcdr3=`cat $entropydir/${brain[$i]}*.all.Hcdr3 | cut -f2`
        Hcdr3_1000=`cat $entropydir/${brain[$i]}*.1000.Hcdr3 | cut -f2`
        Hvj=`cat $entropydir/${brain[$i]}*.all.Hvj | cut -f2`
        Hvj_1000=`cat $entropydir/${brain[$i]}*.1000.Hvj | cut -f2`
        Hv=`cat $entropydir/${brain[$i]}*.all.Hv | cut -f2`
        Hv_1000=`cat $entropydir/${brain[$i]}*.1000.Hv | cut -f2`
        Hj=`cat $entropydir/${brain[$i]}*.all.Hj | cut -f2`
        Hj_1000=`cat $entropydir/${brain[$i]}*.1000.Hj | cut -f2`
	echo -e "${brain[$i]}\tbrain\t${sample_type_br[$i]}\t$reads\t$Hcdr3\t$Hcdr3_1000\t$Hvj\t$Hvj_1000\t$Hv\t$Hv_1000\t$Hj\t$Hj_1000" >> TCR40-43.entropies.tsv
done
