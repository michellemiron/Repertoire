#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <sys/time.h>

using namespace std;

void tokenize(string read, vector<string> &read_vec);

int main(int argc, char** argv){

// read in bam reads
char *sam=argv[1];
string samname=argv[1];
string outdir=argv[2];

ifstream readset;
readset.open(sam);

ofstream merged_readset;
merged_readset.open((outdir+"/"+samname+".merged").c_str());

// read data
string read; // default name to store each line of bam
string head;
int flag; // collect flag associated with a read
string chr; // collect chr information associated with a read
int pos; // starting position on the reference
int len; // end position on the reference
string seq;
string cigar; //collect CIGAR string associated with a read
string quals; // collect quality score associated with a read
string ref;

string prevhead="filler";
while(getline(readset,read)){ // read in one line at a time

vector<string> read_vec;
tokenize(read, read_vec); //tokenize the string
head=read_vec[0];

	if(head==prevhead){
	   flag=atoi(read_vec[1].c_str());
	   //for(int i; i<flag.size();i++) cout<<flag[i]<<endl;
	   chr=read_vec[2];
	   pos=atoi(read_vec[3].c_str());	
	   cigar=read_vec[5];
	   seq=read_vec[9];
	   len=seq.length();
	   quals=read_vec[10];
	   merged_readset<<"\t"<<flag<<" "<<chr<<" "<<pos<<" "<<len<<" "<<seq<<" "<<quals<<" "<<cigar;
	}
	else{
		merged_readset<<endl;
		merged_readset<<head;
		flag=atoi(read_vec[1].c_str());
		chr=read_vec[2];
       		pos=atoi(read_vec[3].c_str());
		seq=read_vec[9];
		len=seq.length();
       		cigar=read_vec[5];
       		quals=read_vec[10];
       		prevhead=head;
		merged_readset<<"\t"<<flag<<" "<<chr<<" "<<pos<<" "<<len<<" "<<seq<<" "<<quals<<" "<<cigar;
	}
}

merged_readset<<endl;
readset.close();
merged_readset.close();


}

/**********************************************/
// other funcs

// tokenize the string
void tokenize(string read, vector<string> &read_vec){

    bool keepEmpty = false;
    
    size_t prev = 0;
    size_t next = 0;	

    while ((next = read.find_first_of("\t", prev)) != string::npos)
    {
        if (keepEmpty || (next - prev != 0))
        {
            read_vec.push_back(read.substr(prev, next - prev));
        }
        prev = next + 1;
    }

    if (prev < read.size())
    {
        read_vec.push_back(read.substr(prev));
    }


}



