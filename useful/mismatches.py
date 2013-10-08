#!/bin/python

import sys
str1=sys.argv[1].strip()
str2=sys.argv[2].strip()

str1=[i for i in str1]
str2=[i for i in str2]


pos=[]
for ind in range(len(str1)):
	if str1[ind]!=str2[ind]:
		pos.append(ind+1)


print pos
print len(pos)
