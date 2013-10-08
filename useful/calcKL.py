import numpy as np
import sys

# Input
	# fbl: frequencies from blood (1st argument)
	# fbr: frequencies from brain (2nd argument)

# read in files
fbl=open(sys.argv[1],'r').readlines()
fbr=open(sys.argv[2],'r').readlines()

# convert to floating point
px=np.array(fbl).astype(float)
qx=np.array(fbr).astype(float)

# compute individual entropies
Hpq=sum(-px*np.log2(qx))
Hp=sum(-px*np.log2(px))

# get KL

print Hpq-Hp

