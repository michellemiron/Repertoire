import numpy as np
import sys

input=sys.stdin.read().strip().split('\n')
input=[i.split('\t') for i in input]
input=zip(*input)

values=np.array(input[1]).astype(float)
px=values/values.sum()
px=px[px.nonzero()]
H=sum(-px*np.log2(px))
print H

