#!/bin/bash

### Comments ###
#
# sed 1d deletes header line
# cut extracts specific columns. First column should be the item, second the value. In the CDR3 case I have an additional 'awk' statement to switch around the columns because cut does not do this, so the sequences comes out after the number.
# merge.awk merges like items in column 1
# entropy.py computes and prints out entropy

#r=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRB/Jenny_out
r=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/fastqcorrected/processedbams/TCRB/Jenny_out
s=$1
i=$r/$s*chain.productive.tsv
#for i in *.productive.tsv;
#do
Hcdr3=`sed 1d $i | cut -f5,2 | awk '$1>1 {print $0}'| awk '{print $2"\t"$1}' |  awk -f merge.awk | python entropy.py`  # for CDR3 amino acid seq  
echo -e "Hcdr3\t$Hcdr3"
Hvj=`sed 1d $i |  cut -f1,2 |awk '$2>1 {print $0}' | awk -f merge.awk | python entropy.py` # for VJ
echo -e "Hvj\t$Hvj"
Hv=`sed 1d $i | cut -f3,2  | awk '$1>1 {print $0}'| awk '{print $2"\t"$1}' | awk -f merge.awk |python entropy.py`  # for CDR3 amino acid seq 
echo -e "Hv\t$Hv"
Hj=`sed 1d $i | cut -f4,2  | awk '$1>1 {print $0}' | awk '{print $2"\t"$1}' | awk -f merge.awk |python entropy.py`  # for CDR3 amino acid seq 
echo -e "Hj\t$Hj"

max=`sed 1d $i | cut -f5,2 | awk '$1>1 {print $0}'|cut -f2 |  sort -u |wc -l`
echo -e "max:\t$max"

#echo -e "$i\t$Hcdr3" >> Hcdr3.txt
#done
