#include <iostream>
#include <fstream>
#include <map>
#include <cstdlib>
#include <string>
#include <cstring>

using namespace std;

string reverse(string seq, int seq_len);
string reverse_complement(string seq, int length);
void read_gene_code(const char *f, map<string, string> &x);
void translate(string ntseq, string &aaseq, map<string,string> codon_table);

map<string, string> codon_table;

int main(int argc, char** argv){

char *nameof_seq = argv[1];
char *nameof_codontable=argv[2];

ifstream stream_seq;
stream_seq.open(nameof_seq);  

string header;
while ( getline(stream_seq,header)){

string seq, qual, revseq, aaseq;
string skip;
string aaseqArray[6];

getline(stream_seq,seq);
getline(stream_seq,skip);
getline(stream_seq,qual);

read_gene_code(nameof_codontable,codon_table);

revseq=reverse_complement(seq,seq.length()); // reverse complement

// read frame 1
translate(seq,aaseq,codon_table);
aaseqArray[0]=aaseq;

// read frame 2
seq=seq.substr(1, seq.length());
translate(seq,aaseq,codon_table);
aaseqArray[1]=aaseq;

// read frame 3
seq=seq.substr(1, seq.length());
translate(seq,aaseq,codon_table); 
aaseqArray[2]=aaseq;

// read frame 1 revcomp
translate(revseq,aaseq,codon_table);
aaseqArray[3]=aaseq;

// read frame 2 revcomp
revseq=revseq.substr(1, revseq.length());
translate(revseq,aaseq,codon_table);
aaseqArray[4]=aaseq;

// read frame 3 revcomp
revseq=revseq.substr(1, revseq.length());
translate(revseq,aaseq,codon_table);
aaseqArray[5]=aaseq;

for (int i=0; i<=5; i++)
	cout<<aaseqArray[i]<<endl;

cout<<endl;







}
stream_seq.close();

}

//////////////////////////////////////////////////
string reverse(string seq, int seq_len){
        string revseq;
        char c;
       for (int i=seq_len-1; i>=0; i--){
           c=seq[i];
        revseq.push_back(c);
        }
        return revseq;
}


string reverse_complement(string seq, int seq_len){
       string revcomp;
       char c;
       for (int i=seq_len-1; i>=0; i--){
           c=seq[i];
           if(c=='A')
                revcomp.push_back('T');
           else if (c=='C')
                revcomp.push_back('G');
           else if (c=='G')
                revcomp.push_back('C');
           else if (c=='T')
                revcomp.push_back('A');
           else
                revcomp.push_back(c);
           }
        return revcomp;
}

// read codonmap.txt (codon \t AA)
void read_gene_code(const char *f, map<string, string> &x)
{
	ifstream is(f);
	string s1, s2;

	x.clear();
	while (1)
	{
		is >> s1;
	if (is.eof()) break;
		is >> s2;
                x[s1] = s2;
        }
}

void translate(string ntseq, string &aaseq, map<string,string> codon_table){

int i=0;

aaseq.clear();
while(ntseq.length()-i>2){
	string codon;
	codon.push_back(char(ntseq[i]));
	codon.push_back(char(ntseq[i+1]));
	codon.push_back(char(ntseq[i+2]));
	i=i+3;
	aaseq+=codon_table[codon];
	}
}
