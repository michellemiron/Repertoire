#!/usr/bin/python

import sys
import re
from copy import deepcopy

# codon dictionary (codon->aa)
codonmap=open('/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/postprocessing/codonmap.txt','r').readlines()
codonmap=[line.strip().split() for line in codonmap]
codondict=(dict(codonmap))

# translate into all allowed reading frames
def translate(seq,codondict):
	allframes=[[],[],[]]
	for rf in range(3):
		pos=deepcopy(rf)
		aatrans=""
		while pos+3<=len(ntseq):
			s=seq[pos:pos+3]
			if re.search('N',s):
				nextaa='-'
			else:
				nextaa=codondict[s]
			aatrans=aatrans+nextaa
			pos=pos+3
		allframes[rf]=aatrans
	return allframes	

# generate sequence output
for line in sys.stdin:
	line=line.strip().split("\t")
	ntseq=line[5].upper()
	aaseq=line[4]
	val=line[1]
	vjpair=line[0]
	vseg=line[2]
	jseg=line[3]
	translations=translate(ntseq,codondict)
	for i,t in enumerate(translations):
		if re.search(aaseq,t):
			m=re.search(aaseq,t)
			aapos=m.span()
			Coffset=1; FGXGoffset=4;
			if (t[aapos[0]-Coffset])=='C': # must have proper C motif
				ntseqCDR3=ntseq[(aapos[0]-Coffset)*3+i:(aapos[1]+FGXGoffset)*3+i] # entire nucleotide sequence from motif start to end
				aaseqCDR3=t[aapos[0]-Coffset:aapos[1]+FGXGoffset] # entire amino acid sequence from motif start to end
				print vjpair+'@'+vseg+'@'+jseg+'@'+aaseqCDR3+'@'+ntseqCDR3+'\t'+val
