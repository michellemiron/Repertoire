#!/bin/python

import sys
#from time import time

def make_kmers(read,k):
	if len(read)>=20:
		for i in range(len(read)-(k-1)):
			yield read[i:i+k] 

#tstart=time()
outfile=sys.argv[1]
f=open(outfile,'a')
for line in sys.stdin:
	line=line.strip()
	allkmers=make_kmers(line,19)
	allkmers=filter(lambda x: str.find(x,"N")==-1,allkmers)	
	for kmer in allkmers:
		f.write(kmer+'\n')

f.close()

#tend=time()
#print 'time_elapsed = ', tend-tstart
#print 'DONE'
