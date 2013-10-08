entropydir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/human/entropy # $1

# human B chain

# blood, blood+brain, GBM, LGG
blood=( TCR27 TCR9 TCR15 TCR13 TCR7 TCR23 TCR25 TCR29 )
sample_type_bl=( blood-healthy1 blood-healthy2 GBM1 GBM2 LGG1 LGG2 LGG3 LGG4 )
brain=( TCR26 TCR14 TCR12 TCR6 TCR22 TCR24 TCR28 )
sample_type_br=( brain-healthy1 GBM1 GBM2 LGG1 LGG2 LGG3 LGG4  )

echo -e "name\ttissue\tsample\tnum_reads\tentropy_CDR3_total\tentropy_CDR3_1000\tnum_VJ\tentropy_VJ_total\tnum_VJ_1000\tentropy_VJ_1000\tnum_V\tentropy_V_total\tnum_V_1000\tentropy_V_1000\tnum_J\tentropy_J_total\tnum_J_1000\tentropy_J_1000">humanBchain.entropies.tsv

n=$[${#blood[@]}-1]
for i in $(seq 0 $n)
do
	reads=`cat  $entropydir/${blood[$i]}*.Bchain*.all.Hcdr3 | cut -f1`
	Hcdr3=`cat  $entropydir/${blood[$i]}*.Bchain*.all.Hcdr3 | cut -f2`
	Hcdr3_1000=`cat $entropydir/${blood[$i]}*.Bchain*.1000.Hcdr3 | cut -f2`
	vjreads=`cat  $entropydir/${blood[$i]}*.Bchain*.all.Hvj | cut -f1`
	Hvj=`cat $entropydir/${blood[$i]}*.Bchain*.all.Hvj | cut -f2`
	vjreads_1000=`cat  $entropydir/${blood[$i]}*.Bchain*.1000.Hvj | cut -f1`
	Hvj_1000=`cat $entropydir/${blood[$i]}*.Bchain*.1000.Hvj | cut -f2`  
        vreads=`cat  $entropydir/${blood[$i]}*.Bchain*.all.Hv | cut -f1`
	Hv=`cat $entropydir/${blood[$i]}*.Bchain*.all.Hv | cut -f2`
	vreads_1000=`cat  $entropydir/${blood[$i]}*.Bchain*.1000.Hv | cut -f1`         
        Hv_1000=`cat $entropydir/${blood[$i]}*.Bchain*.1000.Hv | cut -f2`
	jreads=`cat  $entropydir/${blood[$i]}*.Bchain*.all.Hj | cut -f1`
        Hj=`cat $entropydir/${blood[$i]}*.Bchain*.all.Hj | cut -f2`     
	jreads_1000=`cat $entropydir/${blood[$i]}*.Bchain*.1000.Hj | cut -f1`     
        Hj_1000=`cat $entropydir/${blood[$i]}*.Bchain*.1000.Hj | cut -f2`
	echo -e "${blood[$i]}\tblood\t${sample_type_bl[$i]}\t$reads\t$Hcdr3\t$Hcdr3_1000\t$vjreads\t$Hvj\t$vjreads_1000\t$Hvj_1000\t$vreads\t$Hv\t$vreads_1000\t$Hv_1000\t$jreads\t$Hj\t$jreads_1000\t$Hj_1000" >> humanBchain.entropies.tsv
done 

n=$[${#brain[@]}-1]
for i in $(seq 0 $n)
do
	reads=`cat $entropydir/${brain[$i]}*.Bchain*.all.Hcdr3 | cut -f1`
        Hcdr3=`cat $entropydir/${brain[$i]}*.Bchain*.all.Hcdr3 | cut -f2`
        Hcdr3_1000=`cat $entropydir/${brain[$i]}*.Bchain*.1000.Hcdr3 | cut -f2`
	vjreads=`cat  $entropydir/${brain[$i]}*.Bchain*.all.Hvj | cut -f1`
        Hvj=`cat $entropydir/${brain[$i]}*.Bchain*.all.Hvj | cut -f2`
	vjreads_1000=`cat  $entropydir/${brain[$i]}*.Bchain*.1000.Hvj | cut -f1`
        Hvj_1000=`cat $entropydir/${brain[$i]}*.Bchain*.1000.Hvj | cut -f2`
        vreads=`cat  $entropydir/${brain[$i]}*.Bchain*.all.Hv | cut -f1`
	Hv=`cat $entropydir/${brain[$i]}*.Bchain*.all.Hv | cut -f2`
        vreads_1000=`cat  $entropydir/${brain[$i]}*.Bchain*.1000.Hv | cut -f1`
	Hv_1000=`cat $entropydir/${brain[$i]}*.Bchain*.1000.Hv | cut -f2`
	jreads=`cat  $entropydir/${brain[$i]}*.Bchain*.all.Hj | cut -f1`
        Hj=`cat $entropydir/${brain[$i]}*.Bchain*.all.Hj | cut -f2`
	jreads_1000=`cat  $entropydir/${brain[$i]}*.Bchain*.1000.Hj | cut -f1`
        Hj_1000=`cat $entropydir/${brain[$i]}*.Bchain*.1000.Hj | cut -f2`
	echo -e "${brain[$i]}\tbrain\t${sample_type_br[$i]}\t$reads\t$Hcdr3\t$Hcdr3_1000\t$vjreads\t$Hvj\t$vjreads_1000\t$Hvj_1000\t$vreads\t$Hv\t$vreads_1000\t$Hv_1000\t$jreads\t$Hj\t$jreads_1000\t$Hj_1000" >> humanBchain.entropies.tsv
done


# human A chain

# blood, blood+brain, GBM, LGG
blood=( TCR26 TCR8 TCR14 TCR12 TCR6 TCR22 TCR24 TCR28 ) 
sample_type_bl=( blood-healthy1 blood-healthy2 GBM1 GBM2 LGG1 LGG2 LGG3 LGG4 )
brain=( TCR27 TCR15 TCR13 TCR7 TCR23 TCR25 TCR29 )
sample_type_br=( brain-healthy1 GBM1 GBM2 LGG1 LGG2 LGG3 LGG4  )

echo -e "name\ttissue\tsample\tnum_reads\tentropy_CDR3_total\tentropy_CDR3_1000\tnum_VJ\tentropy_VJ_total\tnum_VJ_1000\tentropy_VJ_1000\tnum_V\tentropy_V_total\tnum_V_1000\tentropy_V_1000\tnum_J\tentropy_J_total\tnum_J_1000\tentropy_J_1000" >humanAchain.entropies.tsv

n=$[${#blood[@]}-1]
for i in $(seq 0 $n)
do
        reads=`cat $entropydir/${blood[$i]}*.Achain*.all.Hcdr3 | cut -f1`
        Hcdr3=`cat $entropydir/${blood[$i]}*.Achain*.all.Hcdr3 | cut -f2`
        Hcdr3_1000=`cat $entropydir/${blood[$i]}*.Achain*.1000.Hcdr3 | cut -f2`
	vjreads=`cat $entropydir/${blood[$i]}*.Achain*.all.Hvj | cut -f1`
        Hvj=`cat $entropydir/${blood[$i]}*.Achain*.all.Hvj | cut -f2`
        vjreads_1000=`cat $entropydir/${blood[$i]}*.Achain*.1000.Hvj | cut -f1`
	Hvj_1000=`cat $entropydir/${blood[$i]}*.Achain*.1000.Hvj | cut -f2`
	vreads=`cat $entropydir/${blood[$i]}*.Achain*.all.Hv | cut -f1`
        Hv=`cat $entropydir/${blood[$i]}*.Achain*.all.Hv | cut -f2`
	vreads_1000=`cat $entropydir/${blood[$i]}*.Achain*.1000.Hv | cut -f1`
        Hv_1000=`cat $entropydir/${blood[$i]}*.Achain*.1000.Hv | cut -f2`
        jreads=`cat $entropydir/${blood[$i]}*.Achain*.all.Hj | cut -f1`
	Hj=`cat $entropydir/${blood[$i]}*.Achain*.all.Hj | cut -f2`
	jreads_1000=`cat $entropydir/${blood[$i]}*.Achain*.1000.Hj | cut -f1`
        Hj_1000=`cat $entropydir/${blood[$i]}*.Achain*.1000.Hj | cut -f2`        
	echo -e "${blood[$i]}\tblood\t${sample_type_bl[$i]}\t$reads\t$Hcdr3\t$Hcdr3_1000\t$vjreads\t$Hvj\t$vjreads_1000\t$Hvj_1000\t$vreads\t$Hv\t$vreads_1000\t$Hv_1000\t$jreads\t$Hj\t$jreads_1000\t$Hj_1000" >> humanAchain.entropies.tsv
done

n=$[${#brain[@]}-1]
for i in $(seq 0 $n)
do
        reads=`cat $entropydir/${brain[$i]}*.Achain*.all.Hcdr3 | cut -f1`
        Hcdr3=`cat $entropydir/${brain[$i]}*.Achain*.all.Hcdr3 | cut -f2`
        Hcdr3_1000=`cat $entropydir/${brain[$i]}*.Achain*.1000.Hcdr3 | cut -f2`
	vjreads=`cat $entropydir/${brain[$i]}*.Achain*.all.Hvj | cut -f1`
        Hvj=`cat $entropydir/${brain[$i]}*.Achain*.all.Hvj | cut -f2`
        vjreads_1000=`cat $entropydir/${brain[$i]}*.Achain*.1000.Hvj | cut -f1`
	Hvj_1000=`cat $entropydir/${brain[$i]}*.Achain*.1000.Hvj | cut -f2`
	vreads=`cat $entropydir/${brain[$i]}*.Achain*.all.Hv | cut -f1`
        Hv=`cat $entropydir/${brain[$i]}*.Achain*.all.Hv | cut -f2`
	vreads_1000=`cat $entropydir/${brain[$i]}*.Achain*.1000.Hv | cut -f1`
        Hv_1000=`cat $entropydir/${brain[$i]}*.Achain*.1000.Hv | cut -f2`
	jreads=`cat $entropydir/${brain[$i]}*.Achain*.all.Hj | cut -f1`
	Hj=`cat $entropydir/${brain[$i]}*.Achain*.all.Hj | cut -f2`
	jreads_1000=`cat $entropydir/${brain[$i]}*.Achain*.1000.Hj | cut -f1`
        Hj_1000=`cat $entropydir/${brain[$i]}*.Achain*.1000.Hj | cut -f2`
        echo -e "${brain[$i]}\tbrain\t${sample_type_br[$i]}\t$reads\t$Hcdr3\t$Hcdr3_1000\t$vjreads\t$Hvj\t$vjreads_1000\t$Hvj_1000\t$vreads\t$Hv\t$vreads_1000\t$Hv_1000\t$jreads\t$Hj\t$jreads_1000\t$Hj_1000" >> humanAchain.entropies.tsv
done
