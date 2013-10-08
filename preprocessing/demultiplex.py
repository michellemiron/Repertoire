#!/usr/bin/python

import sys
from os.path import basename
import re # regex
from copy import deepcopy

processedfile=open(sys.argv[1],'r').readlines()
fname=basename(sys.argv[1])
sequences=open(sys.argv[2],'r').readlines()
barcodes=open(sys.argv[3],'r').readlines()
motif=sys.argv[4]
outdir=sys.argv[5]
barcodelength=6

def make_kmers(seq,k):
	'''
	split sequence into k-mers and store as a generator object
	'''
	for i in range(len(seq)-(k-1)):
		yield seq[i:i+k]

def mutate(code,k):
	'''
	mutate sequence -> allow 1 mismatch at present
	'''
	codelist=[]
	mutlet='ACGT'
	for i in range(k):
		for m in mutlet:
			newcode=code[:i]+m+code[i+1:]
			codelist.append(newcode)
		
	return set(codelist) # does not contain original code

def trim(seq,motif,mutmotifs):
	'''
	trim sequences based on barcode motif
	'''
	newseq=deepcopy(seq); # default output is original sequence
	
	if re.search(motif,seq): # look for expected motif
		newseq=seq.split(motif)[-1] # select sequence after last instance of motif
	else:
		for m in mutmotifs: # list of allowable mutations
			if re.search(m,seq): # regex search for mutated motif
	                	newseq=seq.split(m)[-1] # select sequence after last instance of motif
                        	break;
	return newseq


# construct barcodes and mutate
barcodes=[line.strip().split('\t') for line in barcodes]
barcodedict=dict(barcodes) # dictionary mapping index key to barcode sequence
barcode_mutdict={} # dictionary mapping barcode sequence key to list of all single nt mutations
for ind,code in barcodes:
	codelist=mutate(code,barcodelength) # run mutation code 
	barcode_mutdict[code]=codelist # store mutated sequences

# create and open output files
Fout={}
for key,vals in barcodedict.items():
	exec "%s = open('%s', 'w')" % ('f'+vals, outdir+"/"+fname+"."+vals)

exec "%s = open('%s', 'w')" % ('f_unmatched', outdir+"/"+fname+"."+'unmatched')


# trim sequences according to barcode motif (and all possible mutations)
sequences=[line.strip().split('\t') for line in sequences]
mutmotifs=list(mutate(motif,len(motif))) # mutate barcode motif

tmp=[]
for elem in sequences:
	if(len(elem)<2): # rare anomalous cases where header has no associated sequence (presumably end of read corresponds to J cassette)
		tmp.append(elem)
	else:
		tmp.append([elem[0],trim(elem[1],motif,mutmotifs)])
sequences=tmp;

# use kmer sliding window and match to one of the barcode lists
for i in xrange(len(sequences)):
	ismatched=False;	
	
	# if no sequence present
	if len(sequences[i])<2:
		exec "%s.write(u'''%s''')" % ("f_unmatched",vjline)
		continue;

        ntseq=sequences[i][1]
	vjline=processedfile[i]
	kmerlist=make_kmers(ntseq,barcodelength) # generate windows iteratively
	# look for perfect match	
	while not ismatched:
		try:
			kmer=kmerlist.next()
		except StopIteration:
			kmerlist=make_kmers(ntseq,barcodelength) # remake generator for mutated barcode
			break;
		for key,vals in barcodedict.items():
			if(kmer==vals):
				exec "%s.write(u'''%s''')" % ("f"+vals,vjline) # write output line from pipeline
				ismatched=True
				break;
			
	# look for imperfect (mutated) match
	while not ismatched:
		try:
                        kmer=kmerlist.next()
                except StopIteration:
                        break;
		for key,vals in barcode_mutdict.items():
			if kmer in vals:
				exec "%s.write(u'''%s''')" % ("f"+key,vjline)
				ismatched=True
				break;

	if ismatched is False:
		exec "%s.write(u'''%s''')" % ("f_unmatched",vjline)
		
# close writeable files
for key,vals in barcodedict.items():
        exec "%s.close()" % ("f"+vals)

exec "%s.close()" % ('f_unmatched')
