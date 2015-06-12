__author__ = 'linkuo'
import sys
from string import maketrans
import os
#changer class contains a dict which stores the sequence based on the symbol as key, and a trantab to change "ATCG"to "TAGC"
#has a fullfillsequence method, which take a line of samfile as argument, and change it to bwasw version, using the inside dictionary and the trantab
class changer:
    dict={}
    trantab=0
    def __init__(self, readfile):
        intab = "ATCG"
        outtab="TAGC"
        self.trantab = maketrans(intab, outtab)
        line=readfile.readline()
        symbol=line.split(" ")[0].lstrip("@")
        line=readfile.readline()
        sequence=line.strip()
        self.dict[symbol]=sequence
        for i in range(0,2):
            readfile.readline()
        while 1:
            line=readfile.readline()
            if not line:
                break
            symbol=line.split(' ')[0].lstrip("@")
            line=readfile.readline()
            sequence=line.strip()
            self.dict[symbol]=sequence
            for i in range(0,2):
                readfile.readline()
    def fullfillsequence(self,line):
        line=line.split('\t')
        line[14]=line[14].strip()
        line.append("UNS\n")
        if line[1]=='0' or line[1]=='2048':
            line[1]='0'
            line[5]=line[5].replace('H','S')
            line[9]=self.dict[line[0]]
        if line[1]=='16' or line[1]=='2064':
            line[1]='16'
            line[5]=line[5].replace('H','S')
            line[9]=self.dict[line[0]][::-1].translate(self.trantab)
        return '\t'.join(line)
#argv[1] is the samfile for input. argv[2] is the read file for input. argv[3] is the directory of the sam file
def main():
    outfile=open(sys.argv[3]+"tempfileformodifysam","w")
    readfile=open(sys.argv[2],"r")
    #use the content in readfile to initialize the changer.
    modifier=changer(readfile)
    readfile.close()
    samfile=open(sys.argv[1],"r")
    line=samfile.readline()
    outfile.write(modifier.fullfillsequence(line))
    while 1:
        line=samfile.readline()
        if not line:
            break
        outfile.write(modifier.fullfillsequence(line))
    outfile.close()
    samfile.close()
    os.remove(sys.argv[1])
    os.rename(sys.argv[3]+"tempfileformodifysam",sys.argv[1])


if __name__ == "__main__": main()