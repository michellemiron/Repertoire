#!/usr/bin/env python

# compute a range of diversity measures

import numpy as np
''' 
Available functions:

calc_clonality: Compute Shannon entropy (with Miller-Madow correction for underestimation bias)
calc_simpson: Compute Simpson's index
calc_gini: Compute Gini coefficient
calc_r50: Compute R50 index

'''




def calc_clonality(vals):
	'''
	Compute information entropy (H=-sum(p*log(p))
	Input: 
		ids -> list of identifiers
		vals->list of counts OR frequencies
	'''

	px=vals.astype(float)/vals.sum() # normalize
	px=px[px.nonzero()] # remove 0s
	H=sum(-px*np.log2(px)) # compute entropy
	#N=sum(vals)
	#M=len(vals)
	#H=E_H+(M-1)/(2*N) # Miller-Madow

	Hmax=np.log2(len(px))
	return 1-H/Hmax


def calc_simpson(vals):
	'''
	Simpson's Diversity Index

	'''
	px=vals.astype(float)/vals.sum() # normalize
        px=px[px.nonzero()] # remove 0s
	px=px**2

	return sum(px)


def calc_gini(vals):
	'''
	Gini inequality coefficient

	'''
	vals=sorted(vals) # values
	bins=np.bincount(vals).astype(float) # histogram
	S=sum(bins)
	fy=dict([[i,j/S] for i,j in enumerate(bins) if j>0]) # remove extra vals and map y to pdf f(y)
	y=[i for i,j in enumerate(bins) if j>0] # copy numbers

	wy=[i*fy[i] for i in y] # compute a weight for each each element
	wy.insert(0,0)

	Sn=sum(wy)	
	num=[fy[j]*(sum(wy[0:i-1])+sum(wy[0:i])) for i,j in enumerate(y,start=2)] # f(y)*(Si-1+Si)
	G=1-sum(num)/Sn
	
	return G


def calc_r50(vals):
	'''
	Fraction of counts at which up to 50% of labels are accounted for

	'''
	vals=sorted(vals)
	M=len(vals) # No of labels
	N=sum(vals) # Total counts

	M50=(M/2) # 50% of labels


	keptvals=[i for j,i in enumerate(vals) if (j+1)<=M50]
	keptvals=map(float,keptvals)

	return 1-sum(keptvals)/N



def main():
	import os
	from argparse import ArgumentParser # for specifying input arguments


	def loadfile(fname,h):
		'''
		Loads files. If h flag is True, assumes a header is present and removes it.

		'''
		fin=open(fname,"r")
		
		if h:
			fin.readline()
			data=fin.readlines()
		else:
			data=fin.readlines()

		fin.close()	

		data=[i.strip().split() for i in data]
		data=zip(*data)

		return data[0],np.array(data[1]).astype(int)

	def get_args():
		'''
			Gets input arguments
		'''
			
		parser=ArgumentParser(description="Compute diversity measures.",epilog="Diversity measures can be imported as a module and run from within a python script")
	
		# Optional Arguments
		parser.add_argument("-a", "--header", help="file has a header", action="store_true")
		parser.add_argument("-t","--type", default="clonality",help="diversity measure. Options are: Entropy, Clonality, Simpson, Gini, Hill, R50")
		
		# Positional arguments
		parser.add_argument("filename",help="input filename. Column 1: ids, Column 2: count or frequency")
		
		return parser


	parser=get_args() # get input arguments
	args=parser.parse_args()

	ids,vals=loadfile(args.filename,args.header) # load the file	
	
	assert(len(ids)==len(vals)) # check that ids map to values

	diversity_type=args.type.lower()

	if diversity_type=="clonality": # Shannon Entropy with Miller-Madow bias correction
		D=calc_clonality(vals)
	elif diversity_type=="simpson": # Simpson's Index
		D=calc_simpson(vals)
	elif diversity_type=="gini": # Gini Coefficient
		D=calc_gini(vals)
	elif diversity_type=="r50": # R50
		D=calc_r50(vals)
	else:
		print "Error: unknown diversity measure"

	print D

	
	

if __name__ == '__main__':
        '''
        python diversitymeasures.py -h to learn more

        '''
	main()
