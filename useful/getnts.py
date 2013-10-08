#!/bin/python

import sys

seq=sys.argv[1]
p1=int(sys.argv[2])
p2=int(sys.argv[3])

seq=seq.strip()
s=[i for i in seq]
s=s[p1-1:p2]

print ''.join(s)
