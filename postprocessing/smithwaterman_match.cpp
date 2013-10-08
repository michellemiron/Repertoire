#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cstring>
#include <cmath>
#include <sys/time.h>


using namespace std;

// define functions
double similarity_score(char a,char b);
double find_array_max(double array[],int length);
void insert_at(char arr[], int n, int idx, char val);
void checkfile(int open, char filename[]);
void read_sequence(ifstream &f, string &header, string &seq, string &qual );
string reverse(string seq, int seq_len);
string reverse_complement(string seq, int length);

string datadir="/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/errorcheck2/"; // output directory
int ind;
double mu,delta;

// main script
int main(int argc, char** argv){
clock_t t1, t2;
  t1=clock();
  // read info from arguments
  if(argc!=8){
    cout<<"Give me the proper number of input arguments:"<<endl<<"1 : mu"<<endl;
    cout<<"2 : delta"<<endl<<"3 : filename sequence A"<<endl<<"4 : filename sequence B"<<endl;
    cout<<"5 : maximal length N of sequences"<<endl;
    cout<<"6 : minimum acceptable overlap length"<<endl;
    cout<<"7 : name of output file"<<endl;exit(1);
   }
  mu    = atof(argv[1]); // mismatch penalty
  delta = atof(argv[2]); // gap penalty
  char *nameof_seq_a = argv[3]; //input fastq file 1
  char *nameof_seq_b = argv[4]; //input fastq file 2
   
  int N_max = atoi(argv[5]); // max length
  int cutoff=atoi(argv[6]);  // minimum considered overlap
  char *outprefix=argv[7]; // output prefix

 
  // read the sequences into two vectors:
  ifstream stream_seq_b;          // first define the input-streams for seq_a and seq_b
  stream_seq_b.open(nameof_seq_b);            // the same for seq_b
  checkfile(! stream_seq_b,nameof_seq_b);
  ifstream stream_seq_a;
  stream_seq_a.open(nameof_seq_a);            // open the file for input
  checkfile(! stream_seq_a,nameof_seq_a);     // check, whether the file was opened successfully

// output paths
ofstream matchfile, mismatchfile, rejectfile;
matchfile.open((datadir+outprefix+".matches.fastq").c_str());
mismatchfile.open((datadir+outprefix+".mismatches").c_str());
rejectfile.open((datadir+outprefix+".discard").c_str());

while (stream_seq_a.good()){
  string seq_a,seq_b;
  string header_a,header_b;
  string qual_a,qual_b;

  // read the next sequence in the fastq file
  read_sequence(stream_seq_a, header_a, seq_a, qual_a);
  read_sequence(stream_seq_b, header_b, seq_b, qual_b);

  char *cstr=new char [header_a.length()+1];
  strcpy(cstr, header_a.c_str());
  char * h = strtok(cstr," "); // header id without 01 or 02

  size_t check_a=seq_a.find("REMOVED");
  size_t check_b=seq_b.find("REMOVED");

 if (check_a != string::npos || check_b != string::npos){
	rejectfile<<h<<endl<<"REMOVED"<<endl;        
	continue;
}	

  
  int N_a = seq_a.length();  // get the actual lengths of the sequences
  int N_b = seq_b.length();
  ////////////////////////////////////////////////
  seq_b=reverse_complement(seq_b,N_b); // reverse complement for TCGA sequence
  qual_b=reverse(qual_b,N_b); // reverse quality string
  
  // initialize H
  double H[N_a+1][N_b+1];     
  for(int i=0;i<=N_a;i++){
    for(int j=0;j<=N_b;j++){
      H[i][j]=0.;
    }
  } 

  double temp[4];
  int I_i[N_a+1][N_b+1],I_j[N_a+1][N_b+1];     // Index matrices to remember the 'path' for backtracking
  
  // Smith Waterman

  for(int i=1;i<=N_a;i++){
    for(int j=1;j<=N_b;j++){
      temp[0] = H[i-1][j-1]+similarity_score(seq_a[i-1],seq_b[j-1]); 
      temp[1] = H[i-1][j]-delta;                  
      temp[2] = H[i][j-1]-delta;                 
      temp[3] = 0.;
      H[i][j] = find_array_max(temp,4);
      switch(ind){
      case 0:                                  // score in (i,j) stems from a match/mismatch
   	I_i[i][j] = i-1;
	I_j[i][j] = j-1;
	break;
      case 1:                                  // score in (i,j) stems from a deletion in sequence A
     	I_i[i][j] = i-1;
	I_j[i][j] = j;
	break;
      case 2:                                  // score in (i,j) stems from a deletion in sequence B
      	I_i[i][j] = i;
	I_j[i][j] = j-1;
	break;
      case 3:                                  // (i,j) is the beginning of a subsequence
      	I_i[i][j] = i;
	I_j[i][j] = j;	
	break;
      }
    }
  }

  // search H for the maximal score
  double H_max = 0.;
  int i_max=0,j_max=0;
  for(int i=1;i<=N_a;i++){
    for(int j=1;j<=N_b;j++){
      if(H[i][j]>H_max){
	H_max = H[i][j];
	i_max = i;
	j_max = j;
      }
    }
  }

     // Backtracking from H_max
  int current_i=i_max,current_j=j_max;
  int align_end1=current_i;
  int align_end2=current_j;
  int next_i=I_i[current_i][current_j];
  int next_j=I_j[current_i][current_j];
  int tick=0;
  char consensus_a[N_a+N_b+2],consensus_b[N_a+N_b+2],overlapqual_a[N_a+N_b+2], overlapqual_b[N_a+N_b+2];

  while(((current_i!=next_i) || (current_j!=next_j)) && (next_j!=0) && (next_i!=0)){

    if(next_i==current_i)  {consensus_a[tick] = '-'; overlapqual_a[tick]='#';}                 // deletion in A
    else                   {consensus_a[tick] = seq_a[current_i-1]; overlapqual_a[tick]=qual_a[current_i-1];}  // match/mismatch in A

    if(next_j==current_j)  {consensus_b[tick] = '-'; overlapqual_b[tick]='#';}                // deletion in B
    else                   {consensus_b[tick] = seq_b[current_j-1]; overlapqual_b[tick]=qual_b[current_i-1];}  // match/mismatch in B

    current_i = next_i;
    current_j = next_j;
    next_i = I_i[current_i][current_j];
    next_j = I_j[current_i][current_j];
    tick++;
    }
 
  string ca,cb,qa,qb;
  for(int i=tick-1;i>=0;i--) {ca.push_back(consensus_a[i]); qa.push_back(overlapqual_a[i]);}
  for(int j=tick-1;j>=0;j--) {cb.push_back(consensus_b[j]); qb.push_back(overlapqual_b[j]);}

//  cout<<header_a<<"\t"<<header_b;


// remove N nucleotides
string tmp;
for(int i=0;i<ca.length();i++){
	if(ca[i]!='N')
	tmp.push_back(ca[i]);
}

// remove any overlaps with length below threshhold value (after removal of Ns)
if(tmp.length()<cutoff)
      rejectfile<<h<<endl<<ca<<"\t"<<cb<<endl;	
else if(ca==cb)
	matchfile<<h<<endl<<ca<<endl<<'+'<<endl<<qa<<endl;
else if(ca!=cb)
	mismatchfile<<h<<endl<<ca<<"\t"<<cb<<endl;
else
	cout<<"Weirdness in: "<<h<<endl;
//else if (ca==cb){ // find matches 
//	c_len=ca.length();
//	}

  
/*
  cout <<endl;
  cout<<align_end1-1<<'\t';
  cout<<align_end1-tick<<'\t';
  cout<< align_end2-1<<'\t';
  cout<<align_end2-tick;
  cout<<endl<<endl;
  */


}

stream_seq_a.close();
stream_seq_b.close();
matchfile.close();
mismatchfile.close();
rejectfile.close();

t2=clock();
float t=float(t2-t1)/CLOCKS_PER_SEC;
cout<<t<<endl;
} // END of main



/////////////////////////////////////////////////////////////////////////////
// auxiliary functions used by main:
/////////////////////////////////////////////////////////////////////////////


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

/////////////////////////////////////////////////////////////////////////////

void checkfile(int open, char filename[]){

  if (open){cout << "Error: Can't open the file "<<filename<<endl;exit(1);}
  //else cout<<"Opened file "<<filename<<endl;
}

/////////////////////////////////////////////////////////////////////////////

double similarity_score(char a,char b){

  double result;
  if(a==b){
      result=1.;
    }
  else{
      result=-mu;
    }
  return result;
}

/////////////////////////////////////////////////////////////////////////////

double find_array_max(double array[],int length){

  double max = array[0];            // start with max = first element
  ind = 0;

  for(int i = 1; i<length; i++){
      if(array[i] > max){
	max = array[i];
	ind = i; 
      }
  }
  return max;                    // return highest value in array
}

/////////////////////////////////////////////////////////////////////////////
void read_sequence(ifstream &f, string &header, string &seq, string &qual )
{
char h[5000];
char seqline[5000]; 
char skipline[10];
char q[5000];
  
  f.getline(h,5000);
  f.getline(seqline,5000);
  f.getline(skipline,10);
  f.getline(q,5000); 

  header=string(h);
  qual=string(q);
      for(int i = 0; seqline[i] != 0; ++i)
	{
	  int c = toupper(seqline[i]);
	  seq.push_back(char(c));
	}
}
