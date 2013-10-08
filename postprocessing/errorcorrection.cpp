#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <cstring>
#include <cmath>
#include <sys/time.h>
#include <smithwaterman.h>

using namespace std;

void tokenize(string read, vector<string> &read_vec, string delim);
string searchfasta(string chr, int pos1,int pos2, string fastapath);
//string reverse(string seq, int seq_len);
string reverse_complement(string seq, int length);
void parse_cigar(string cigar, vector<int> &cigar_num, vector<char> &cigar_let );
void get_matchseq(string seq,string qual, string chr, int pos, vector<int> cigar_num, vector<char> cigar_let,string fastapath, string &refseq, int &lclip);

int main(int argc, char** argv){
// read in bam reads
char *sam1=argv[1];
char *sam2=argv[2];
string fastapath=argv[3];
string outpath=argv[4];
string outname=argv[5];

ifstream readset1, readset2;
readset1.open(sam1);
readset2.open(sam2);

// output paths
ofstream output;
output.open((outpath+"/"+outname+".fastq").c_str());


// read data
string read1; // default name to store each line of bam
string read2;
string head;

getline(readset1,read1); // remove blank first line
getline(readset2,read2);
//int c=0;
while(getline(readset1,read1)){ // read in one line at a time
//c++; cout<<c<<endl;
vector<string> read1_vec;
vector<string> read2_vec;
getline(readset2,read2);

tokenize(read1,read1_vec,"\t");
tokenize(read2,read2_vec,"\t");
head=read1_vec[0];

/* cases
1) Two entries for each read
2) Two entries for read 1, three or four entries for read 2
3) Three or four entries for read 1, two entries for read 2
4) Three or four entries for both reads

Note: entry 1 is always the read id
*/

int rz1=read1_vec.size();
int rz2=read2_vec.size();
string fixedread; // overlapping region post correction 
string fixedqual; // quality score for overlapping region post correction
string finalread="N", finalqual="#"; // read and quality of completely merged reads (overlap and nonoverlap combined). Default is N,# to indicate that a suitable output read could not be produced from the two input reads.
string errorcase="other"; // keep track of the type of correction 0/0,0/1,1/1,2/0,2/1,2/2 to indicate number of usable mapped regions in the two reads.
if(rz1==2 && rz2==2){ // one site mapped for both reads. merge where possible to get longer read that can map in two places.
	errorcase="1/1dump";
	string r1=read1_vec[1];
	string r2=read2_vec[1];
	vector<string> map1_vec;
	vector<string> map2_vec;

	tokenize(r1, map1_vec, " ");
	tokenize(r2, map2_vec, " ");

	int flag1=atoi(map1_vec[0].c_str()), flag2=atoi(map2_vec[0].c_str());
	string chr1=map1_vec[1], chr2=map2_vec[1];
	int pos1=atoi(map1_vec[2].c_str()), pos2=atoi(map2_vec[2].c_str());
	int len1=atoi(map1_vec[3].c_str()), len2=atoi(map2_vec[3].c_str());	
	string seq1=map1_vec[4]; string seq2=map2_vec[4];
	string qual1=map1_vec[5]; string qual2=map2_vec[5];
	
	if (flag1==4 and flag2==4) // discard unmapped reads
		errorcase="0/0";

	else if ((flag1==4 and flag2!=4) or (flag1!=4 and flag2==4)) // discard case where there is only one mapped region
		errorcase="1/0";
	//merge regions with only one mapped segment each
	else if(chr1==chr2 and abs(pos1-pos2)>=250){ // check that mapped segments are different and on same chromosome
		errorcase="1/1fix";
		string ca, cb; // overlapping strings after running smith waterman
		string qa, qb; // quality scores corresponding to overlaps
		int starta, startb, olen; // location on strings where overlap occurs
	
		runsw(1, 4, seq1, seq2, qual1, qual2, 250, ca, cb, qa, qb, starta, startb,olen); // smith waterman alignment
		// ca,cb -- consensus sequence
		// qa,qb -- quality for consensus sequence
		// starta,startb -- starting position of overlap (at the end of each read. right most coordinate)
		// olen -- length of overlap (i.e. length of consensus)
		
		// construct fixed read
	        for(int i=0; i<olen; i++){
			if(ca[i] == cb[i])
			{
			   fixedread.push_back(ca[i]);
			   fixedqual.push_back(qa[i]); 
			}
			else
			{
			  if(int(qa[i])-int(qb[i])>10){fixedread.push_back(ca[i]); fixedqual.push_back(qa[i]);}
			  else if (int(qb[i])-int(qa[i])>10) {fixedread.push_back(cb[i]); fixedqual.push_back(qb[i]);}
			  else {fixedread.push_back(ca[i]); fixedqual.push_back(qa[i]);}
			}
		}
		
		// get nonoverlapping parts and attach to overlap
		string seqstart, seqend, qualstart, qualend;
		int cutoff_left, cutoff_right, lastind;// start and end for fixed read using nonoverlapping regions
		if(startb-olen<starta-olen) {seqstart=seq1; qualstart=qual1; cutoff_left=starta-olen;} else {seqstart=seq2; qualstart=qual2; cutoff_left=startb-olen;}
		if((len1-starta)<(len2-startb)) {seqend=seq2; qualend=qual2; cutoff_right=startb; lastind=len2;}  else {seqend=seq1; qualend=qual1; cutoff_right=starta; lastind=len1;}

		// merge
		string frontpart, frontqual;
		for(int i=0;i<=cutoff_left;i++){
			frontpart.push_back(seqstart[i]);
			frontqual.push_back(qualstart[i]);
		}
		string backpart, backqual;
		for(int i=cutoff_right+1;i<=lastind;i++){
                        backpart.push_back(seqend[i]);
			backqual.push_back(qualend[i]);
                }

		finalread.append(frontpart); finalread.append(fixedread); finalread.append(backpart);
		finalqual.append(frontqual); finalqual.append(fixedqual); finalqual.append(backqual);
		//cout<<abs(starta-olen-(startb-olen))<<endl;
		//cout<<seq1<<"\t"<<seq1.length()<<endl<<seq2<<"\t"<<seq2.length()<<endl<<finalread<<"\t"<<finalread.length()<<endl<<endl;
		//cout<<finalread<<endl<<finalqual<<endl<<endl;
		//cout<<head<<"\t"<<finalread<<"\t"<<pos1<<"\t"<<pos2<<endl<<endl;
	}} // end 1-1 comparison

else if((rz1==3 or rz1==4) && rz2==2){ // one read mapped twice or more, the other only once. use the singly mapped read to correct the double mapped.
	errorcase="2/1dump";
 	string r1a=read1_vec[1];
	string r1b=read1_vec[2];
        string r2=read2_vec[1];

        vector<string> map1a_vec;
	vector<string> map1b_vec;
        vector<string> map2_vec;
        tokenize(r1a, map1a_vec, " ");
	tokenize(r1b, map1b_vec, " ");
        tokenize(r2, map2_vec, " ");

        int flag1a=atoi(map1a_vec[0].c_str()), flag1b=atoi(map1a_vec[0].c_str()), flag2=atoi(map2_vec[0].c_str());
        string chr1a=map1a_vec[1], chr1b=map1b_vec[1], chr2=map2_vec[1];
        int pos1a=atoi(map1a_vec[2].c_str()), pos1b=atoi(map1b_vec[2].c_str()), pos2=atoi(map2_vec[2].c_str());
        //int len1=atoi(map1_vec[3].c_str()), len2=atoi(map2_vec[3].c_str());
        string seq1=map1a_vec[4], seq2=map2_vec[4];
        string qual1=map1a_vec[5], qual2=map2_vec[5];
	string cigar1a=map1a_vec[6], cigar1b=map1b_vec[6], cigar2=map2_vec[6];
	
 	if(rz1==4) { // mapped to 3 places. check which mapped regions belong together and make them a and b.
		vector<string> map1c_vec;
		string r1c=read1_vec[3]; tokenize(r1c, map1c_vec, " ");
                if (chr1a==map1c_vec[1] and chr1a!=chr1b) {flag1b=atoi(map1c_vec[0].c_str()); chr1b=map1c_vec[1]; pos1b=atoi(map1c_vec[2].c_str()); cigar1b=map1c_vec[6];}
                else if (chr1b==map1c_vec[1] and chr1a!=chr1b){flag1a=atoi(map1c_vec[0].c_str()); chr1a=map1c_vec[1]; pos1a=atoi(map1c_vec[2].c_str()); cigar1a=map1c_vec[6];}
		}
	if(flag2==4){
		errorcase="2/0";
		finalread=seq1;
		finalqual=qual1;
	}
	else {
		if(chr1a==chr1b and chr1a==chr2){ // makes sure the reads and mapping are reliable
		// use the singly mapped read to correct the double mapped read
		string ca,cb,qa,qb,mcigar; int starta,startb,olen,mpos,mlen,mflag; int lclip=0;
                vector<char> ciga_let, cigb_let;
                vector<int> ciga_num, cigb_num;		
		int matchpos, refpos; string refseq=""; // starting read position of interest, reference position relative to overlap, and reference sequence
		
		if(abs(pos1a-pos2) < abs(pos1b-pos2)) {mpos=pos1a; mflag=flag1a; mcigar=cigar1a; mpos=pos1a;}
                else {mpos=pos1b; mflag=flag1b;  mcigar=cigar1b; mpos=pos1b;}
		//cout<<head<<"\t"<<mcigar<<"\t";
		if (mflag!=flag2  and abs(pos2-mpos)<=200){
		errorcase="2/1fix";
		parse_cigar(mcigar,ciga_num,ciga_let);  parse_cigar(cigar2,cigb_num,cigb_let);  // parse cigar strings
		get_matchseq(seq1,qual1,chr1a,mpos,ciga_num, ciga_let,fastapath,refseq,lclip); // get reference sequence. lclip is the length of the clipping on the left end of the sequence (to know where the match starts).
		runsw(1, 4, seq1, seq2, qual1, qual2, 250, ca, cb, qa, qb, starta, startb,olen); // get overlap
		
		// check if we can get additional nt from the ends of the string
		 string sloverhang=""; string  qloverhang=""; string sroverhang=""; string qroverhang="";
                if(startb>starta){
                        for(int i=0; i<(startb-starta);i++)
                        {
                                sloverhang.push_back(seq2[i]); qloverhang.push_back(qual2[i]);
                        }
                }
                if((seq2.length()-startb)>(seq1.length()-starta))
                {
                        for(int i=seq1.length();i<seq2.length();i++)
                        {
                                sroverhang.push_back(seq2[i]); qroverhang.push_back(qual2[i]);
                        }
                }
		
		string frontpart, frontqual, backpart,backqual; // nonoverlapping parts
		for(int i=0;i<=starta-olen;i++){
			frontpart.push_back(seq1[i]);
			frontqual.push_back(qual1[i]);
		}
		for(int i=starta+1;i<seq1.length();i++){
			backpart.push_back(seq1[i]);
			backqual.push_back(qual1[i]);
		}

		// overlap
		int j=0; // keep track of reference nts
		int k=startb-olen+1; // overlapping indices for sequence 2
		string overlappart, overlapqual;
		for(int i=starta-olen+1; i<=starta; i++){
			if(seq1[i]==seq2[k] ) {overlappart.push_back(seq1[i]); overlapqual.push_back(qual1[i]);}
			else
			{
				if(i<lclip or i>=(lclip+refseq.length()))
				{ // where there is no reference, go by quality
                          		if(int(qual1[i])-int(qual2[k])>10){overlappart.push_back(seq1[i]); overlapqual.push_back(qual1[i]);}
                          		else if (int(qual2[k])-int(qual1[i])>10) {overlappart.push_back(seq2[k]); overlapqual.push_back(qual2[k]);}
                          		else {overlappart.push_back(seq1[i]); overlapqual.push_back(qual1[i]);}				
				}
				else { // use reference where possible
					j=i-lclip; // place in reference
					if (seq1[i]==refseq[j]) {overlappart.push_back(seq1[i]); overlapqual.push_back(qual1[i]);} // use matching reference
					else if (seq2[k]==refseq[j]) {overlappart.push_back(seq2[k]); overlapqual.push_back(qual2[k]);} // use matching reference
					else {if(seq1[i]=='N'){overlappart.push_back(seq2[k]); overlapqual.push_back(qual2[k]);} else {overlappart.push_back(seq1[i]); overlapqual.push_back(qual1[i]);}} // where neither matches reference use sequence 1, unless sequence 1 is N
				}
			}
			k++;
		}
		finalread.append(sloverhang); finalread.append(frontpart); finalread.append(overlappart); finalread.append(backpart); finalread.append(sroverhang);
		finalqual.append(qloverhang); finalqual.append(frontqual); finalqual.append(overlapqual); finalqual.append(backqual); finalqual.append(qroverhang);
}}}} 


else if(rz1==2 && (rz2==3 or rz2==4)){ // one read mapped twice or more, the other only once. use the singly mapped read to correct the double mapped (reverse case)
	errorcase="2/1dump";
	string r1a=read2_vec[1]; // same as previous else if, but flip the reads and do the same analysis
	string r1b=read2_vec[2];
        string r2=read1_vec[1];

        vector<string> map1a_vec;
	vector<string> map1b_vec;
        vector<string> map2_vec;
        tokenize(r1a, map1a_vec, " ");
	tokenize(r1b, map1b_vec, " ");
        tokenize(r2, map2_vec, " ");

        int flag1a=atoi(map1a_vec[0].c_str()), flag1b=atoi(map1a_vec[0].c_str()), flag2=atoi(map2_vec[0].c_str());
        string chr1a=map1a_vec[1], chr1b=map1b_vec[1], chr2=map2_vec[1];
        int pos1a=atoi(map1a_vec[2].c_str()), pos1b=atoi(map1b_vec[2].c_str()), pos2=atoi(map2_vec[2].c_str());
        //int len1=atoi(map1_vec[3].c_str()), len2=atoi(map2_vec[3].c_str());
        string seq1=map1a_vec[4], seq2=map2_vec[4];
        string qual1=map1a_vec[5], qual2=map2_vec[5];
	string cigar1a=map1a_vec[6], cigar1b=map1b_vec[6], cigar2=map2_vec[6];

       if(rz2==4) { // mapped to 3 places. check which mapped regions belong together and make them a and b.
                vector<string> map1c_vec; string r1c=read2_vec[3];  tokenize(r1c, map1c_vec, " ");
                if (chr1a==map1c_vec[1] and chr1a!=chr1b) {flag1b=atoi(map1c_vec[0].c_str()); chr1b=map1c_vec[1]; pos1b=atoi(map1c_vec[2].c_str()); cigar1b=map1c_vec[6];}
                else if (chr1b==map1c_vec[1] and chr1a!=chr1b){flag1a=atoi(map1c_vec[0].c_str()); chr1a=map1c_vec[1]; pos1a=atoi(map1c_vec[2].c_str()); cigar1a=map1c_vec[6];}
	}
	if(flag2==4){
		finalread=seq1;
		finalqual=qual1;
		errorcase="2/0";
	}
	else {
		if(chr1a==chr1b and chr1a==chr2){ // makes sure the reads and mapping are reliable
		// use the singly mapped read to correct the double mapped read
		string ca,cb,qa,qb,mcigar; int starta,startb,olen,mpos,mlen,mflag; int lclip=0;
                vector<char> ciga_let, cigb_let;
                vector<int> ciga_num, cigb_num;		
		int matchpos, refpos; string refseq=""; // starting read position of interest, reference position relative to overlap, and reference sequence
		
		if(abs(pos1a-pos2) < abs(pos1b-pos2)) {mpos=pos1a; mflag=flag1a; mcigar=cigar1a; mpos=pos1a;}
                else {mpos=pos1b; mflag=flag1b; mcigar=cigar1b; mpos=pos1b;}
		//cout<<head<<"\t"<<mcigar<<"\t";
		
		if (mflag!=flag2  and abs(pos2-mpos)<=200){
		errorcase="2/1fix";
		parse_cigar(mcigar,ciga_num,ciga_let);  parse_cigar(cigar2,cigb_num,cigb_let);  // parse cigar strings
		get_matchseq(seq1,qual1,chr1a,mpos,ciga_num, ciga_let,fastapath,refseq,lclip); // get reference sequence
		runsw(1, 4, seq1, seq2, qual1, qual2, 250, ca, cb, qa, qb, starta, startb,olen); // get overlap	
		// check if we can get additional nt from the ends of the string
		 string sloverhang=""; string  qloverhang=""; string sroverhang=""; string qroverhang="";
                if(startb>starta){
                        for(int i=0; i<(startb-starta);i++)
                        {
                                sloverhang.push_back(seq2[i]); qloverhang.push_back(qual2[i]);
                        }
                }
                if((seq2.length()-startb)>(seq1.length()-starta))
                {
                        for(int i=seq1.length();i<seq2.length();i++)
                        {
                                sroverhang.push_back(seq2[i]); qroverhang.push_back(qual2[i]);
                        }
                }
		
		string frontpart, frontqual, backpart,backqual; // nonoverlapping parts
		for(int i=0;i<=starta-olen;i++){
			frontpart.push_back(seq1[i]);
			frontqual.push_back(qual1[i]);
		}
		for(int i=starta+1;i<seq1.length();i++){
			backpart.push_back(seq1[i]);
			backqual.push_back(qual1[i]);
		}

		// overlap
		int j=0; // keep track of reference nts
		int k=startb-olen+1; // overlapping indices for sequence 2
		string overlappart, overlapqual;
		for(int i=starta-olen+1; i<=starta; i++){
			if(seq1[i]==seq2[k] ) {overlappart.push_back(seq1[i]); overlapqual.push_back(qual1[i]);}
			else
			{
				if(i<lclip or i>=(lclip+refseq.length()))
				{ // where there is no reference, go by quality
                          		if(int(qual1[i])-int(qual2[k])>10){overlappart.push_back(seq1[i]); overlapqual.push_back(qual1[i]);}
                          		else if (int(qual2[k])-int(qual1[i])>10) {overlappart.push_back(seq2[k]); overlapqual.push_back(qual2[k]);}
                          		else {overlappart.push_back(seq1[i]); overlapqual.push_back(qual1[i]);}				
				}
				else { // use reference where possible
					j=i-lclip; // place in reference
					if (seq1[i]==refseq[j]) {overlappart.push_back(seq1[i]); overlapqual.push_back(qual1[i]);} // use matching reference
					else if (seq2[k]==refseq[j]) {overlappart.push_back(seq2[k]); overlapqual.push_back(qual2[k]);} // use matching reference
					else {if(seq1[i]=='N'){overlappart.push_back(seq2[k]); overlapqual.push_back(qual2[k]);} else {overlappart.push_back(seq1[i]); overlapqual.push_back(qual1[i]);}} // where neither matches reference use sequence 1, unless sequence 1 is N
				}
			}
			k++;
		}
		finalread.append(sloverhang); finalread.append(frontpart); finalread.append(overlappart); finalread.append(backpart); finalread.append(sroverhang);
		finalqual.append(qloverhang); finalqual.append(frontqual); finalqual.append(overlapqual); finalqual.append(backqual); finalqual.append(qroverhang);
		//cout<<seq1<<endl<<seq2<<endl<<finalread<<endl<<endl;	
}}}}

else if((rz1==3 or rz1==4) and (rz2==3 or rz2==4)){ // both are mapped
	errorcase="2/2dump";
        string r1a=read1_vec[1]; // same as previous else if, but flip the reads and do the same analysis
        string r1b=read1_vec[2];
        string r2a=read2_vec[1];
	string r2b=read2_vec[2];
        vector<string> map1a_vec;
        vector<string> map1b_vec;
	vector<string> map1c_vec;
        vector<string> map2a_vec;
	vector<string> map2b_vec;
	vector<string> map2c_vec;
        tokenize(r1a, map1a_vec, " ");
        tokenize(r1b, map1b_vec, " ");
        tokenize(r2a, map2a_vec, " ");
	tokenize(r2b, map2b_vec, " ");
        int flag1a=atoi(map1a_vec[0].c_str()), flag1b=atoi(map1b_vec[0].c_str()), flag2a=atoi(map2a_vec[0].c_str()), flag2b=atoi(map2b_vec[0].c_str());
        string chr1a=map1a_vec[1], chr1b=map1b_vec[1], chr2a=map2a_vec[1], chr2b=map2b_vec[1];
        int pos1a=atoi(map1a_vec[2].c_str()), pos1b=atoi(map1b_vec[2].c_str()), pos2a=atoi(map2a_vec[2].c_str()), pos2b=atoi(map2b_vec[2].c_str());
        //int len1=atoi(map1_vec[3].c_str()), len2=atoi(map2_vec[3].c_str());
        string seq1=map1a_vec[4], seq2=map2a_vec[4];
        string qual1=map1a_vec[5], qual2=map2a_vec[5];
        string cigar1a=map1a_vec[6], cigar1b=map1b_vec[6], cigar2a=map2a_vec[6], cigar2b=map2b_vec[6];

        if(rz1==4) {  // mapped to 3 places. check which mapped regions belong together and make them a and b.
		string r1c=read1_vec[3]; tokenize(r1c, map1c_vec, " ");
                if (chr1a==map1c_vec[1] and chr1a!=chr1b) {flag1b=atoi(map1c_vec[0].c_str()); chr1b=map1c_vec[1]; pos1b=atoi(map1c_vec[2].c_str()); cigar1b=map1c_vec[6];}
                else if (chr1b==map1c_vec[1] and chr1a!=chr1b){flag1a=atoi(map1c_vec[0].c_str()); chr1a=map1c_vec[1]; pos1a=atoi(map1c_vec[2].c_str()); cigar1a=map1c_vec[6];}
        }
        if(rz2==4){ // mapped to 3 places. check which mapped regions belong together and make them a and b.
		string r2c=read2_vec[3]; tokenize(r2c, map2c_vec, " ");
                if (chr2a==map2c_vec[1] and chr2a!=chr2b) {flag2b=atoi(map2c_vec[0].c_str()); chr2b=map2c_vec[1]; pos2b=atoi(map2c_vec[2].c_str()); cigar2b=map2c_vec[6];}
                else if (chr2b==map2c_vec[1] and chr2a!=chr2b){flag2a=atoi(map2c_vec[0].c_str()); chr2a=map2c_vec[1]; pos2a=atoi(map2c_vec[2].c_str()); cigar2a=map2c_vec[6];}
	}

	string ca,cb,qa,qb,mcigar; int starta,startb,olen,mpos,mlen,mflag; int lclip=0;
	if (chr1a==chr1b and chr1a==chr2a and chr1a==chr2b){
		runsw(1, 4, seq1, seq2, qual1, qual2, 250, ca, cb, qa, qb, starta, startb,olen); // get smith waterman alignment
                        if(ca==cb){errorcase="2/2match"; finalread=seq1; finalqual=qual1;} // alignment is exact with no mismatches
			else{  // first determine which mapped regions are on the left end and which on the right for each read so that the right regions are compared and pthe final read is correctly assembled	
				int leftpos1,rightpos1,leftpos2,rightpos2,lclip1left=0,lclip2left=0, lclip1right=0, lclip2right=0, leftflag1, rightflag1, leftflag2,rightflag2;
				string refseq1left,refseq1right,refseq2left,refseq2right;
				vector<char> ciga1_let, cigb1_let,ciga2_let, cigb2_let; vector<int> ciga1_num, cigb1_num, ciga2_num, cigb2_num;
				string leftcigar1,leftcigar2,rightcigar1,rightcigar2;  // cigar strings			
				string refseq_left_use, refseq_right_use,sequse, qualuse,  seqcompare, qualcompare; // one string to keep, one to use for comparison
				int start1, start2, lclipleft_use, lclipright_use; // clipping and overlap starting positions (once main string is identified)
				
				if(pos1a<pos1b){leftpos1=pos1a; rightpos1=pos1b; leftcigar1=cigar1a; rightcigar1=cigar1b; leftflag1=flag1a; rightflag1=flag1b;}
				 else{leftpos1=pos1b; rightpos1=pos1a; leftcigar1=cigar1b; rightcigar1=cigar1a; leftflag1=flag1b; rightflag1=flag1a;} //determine which pos is first
				 if(pos2a<pos2b){leftpos2=pos2a; rightpos2=pos2b; leftcigar2=cigar2a; rightcigar2=cigar2b; leftflag2=flag2a; rightflag2=flag2b;} 
				else{leftpos2=pos2b; rightpos2=pos2a; leftcigar2=cigar2b; rightcigar2=cigar2a; leftflag2=flag2b; rightflag2=flag2a;}
				
				// left case
				if(leftflag1!=leftflag2 and rightflag1!=rightflag2){
				errorcase="2/2fix";
				parse_cigar(leftcigar1,ciga1_num,ciga1_let);  parse_cigar(leftcigar2,cigb1_num,cigb1_let);  // parse cigar strings	
				get_matchseq(seq1,qual1,chr1a,leftpos1,ciga1_num, ciga1_let,fastapath,refseq1left,lclip1left); // get reference sequence	
				get_matchseq(seq2,qual2,chr1a,leftpos2,cigb1_num, cigb1_let,fastapath,refseq2left,lclip2left); // get reference sequence
				//cout<<ciga_num[1]<<endl<<cigb_num[1]<<endl;
				//right case 
				parse_cigar(rightcigar1,ciga2_num,ciga2_let);  parse_cigar(rightcigar2,cigb2_num,cigb2_let);  // parse cigar strings  
                                get_matchseq(seq1,qual1,chr1a,rightpos1,ciga2_num, ciga2_let,fastapath,refseq1right,lclip1right); // get reference sequence  
                                get_matchseq(seq2,qual2,chr1a,rightpos2,cigb2_num, cigb2_let,fastapath,refseq2right,lclip2right); // get reference sequence
				
				//cout<<leftcigar1<<endl;//<<leftcigar2<<endl<<endl;
				//cout<<rightcigar1<<endl<<endl; //<<rightcigar2<<endl<<endl;
				//cout<<refseq2left<<endl<<refseq2right<<endl<<endl;
				//cout<<chr1a<<endl<<leftcigar1<<endl<<leftpos1<<endl<<rightcigar1<<endl<<rightpos1<<endl<<refseq1left<<endl<<refseq1right<<endl<<endl;
				if((refseq1left.length()+refseq1right.length())>(refseq2left.length()+refseq2right.length())){	// use as main sequence the one with the most matches
					refseq_left_use=refseq1left; refseq_right_use=refseq1right; sequse=seq1; seqcompare=seq2; qualuse=qual1; qualcompare=qual2; lclipleft_use=lclip1left; lclipright_use=lclip1right; start1=starta; start2=startb;
				}
				else{refseq_left_use=refseq2left; refseq_right_use=refseq2right; sequse=seq2; qualuse=qual2; seqcompare=seq1; qualcompare=qual1; lclipleft_use=lclip2left; lclipright_use=lclip2right; start1=startb; start2=starta;} 
				
			string frontpart, frontqual, overlap, backpart,backqual; // nonoverlapping parts
			for(int i=0;i<=start1-olen;i++){
				frontpart.push_back(sequse[i]);
				frontqual.push_back(qualuse[i]);
			}
			for(int i=start1+1;i<sequse.length();i++){
				backpart.push_back(sequse[i]);
				backqual.push_back(qualuse[i]);
			}

		// overlap
		int j=0; // keep track of reference nts
		int k=start2-olen+1; // overlapping indices
		string overlappart, overlapqual;
		for(int i=start1-olen+1; i<=start1; i++)
		{
			if(sequse[i]==seqcompare[k] ) {overlappart.push_back(sequse[i]); overlapqual.push_back(qualuse[i]);}
			else
			{       // nonreference
				if(i<lclipleft_use or (i>=(lclipleft_use+refseq_left_use.length()) and i<lclipright_use) or (i>=(lclipright_use+refseq_right_use.length())))
				{ // where there is no reference, go by quality
                          		if(int(qualuse[i])-int(qualcompare[k])>10){overlappart.push_back(sequse[i]); overlapqual.push_back(qualuse[i]);}
                          		else if (int(qualcompare[k])-int(qualuse[i])>10) {overlappart.push_back(seqcompare[k]); overlapqual.push_back(qualcompare[k]);}
                          		else {overlappart.push_back(sequse[i]); overlapqual.push_back(qualuse[i]);}				
				}
				else if (i>=lclipleft_use and i<(lclipleft_use+refseq_left_use.length())) { // use reference for match 1
					j=i-lclipleft_use; // place in reference
					if (sequse[i]==refseq_left_use[j]) {overlappart.push_back(sequse[i]); overlapqual.push_back(qualuse[i]);} // use matching reference
					else if (seqcompare[k]==refseq_left_use[j]) {overlappart.push_back(seqcompare[k]); overlapqual.push_back(qualcompare[k]);} // use matching reference
					else {if(sequse[i]=='N'){overlappart.push_back(seqcompare[k]); overlapqual.push_back(qualcompare[k]);} else {overlappart.push_back(sequse[i]); overlapqual.push_back(qualuse[i]);}} // where neither matches reference use sequence 1, unless sequence 1 is N
				}
				else if (i>=lclipright_use and i<(lclipright_use+refseq_right_use.length())) { // use reference for match2
                                        j=i-lclipright_use; // place in reference
                                        if (sequse[i]==refseq_right_use[j]) {overlappart.push_back(sequse[i]); overlapqual.push_back(qualuse[i]);} // use matching reference
                                        else if (seqcompare[k]==refseq_right_use[j]) {overlappart.push_back(seqcompare[k]); overlapqual.push_back(qualcompare[k]);} // use matching reference
                                        else {if(sequse[i]=='N'){overlappart.push_back(seqcompare[k]); overlapqual.push_back(qualcompare[k]);} else {overlappart.push_back(sequse
[i]); overlapqual.push_back(qualuse[i]);}} // where neither matches reference use sequence 1, unless sequence 1 is N
                                }
			}
			k++;
		}
		 finalread.append(frontpart); finalread.append(overlappart); finalread.append(backpart);
		 finalqual.append(frontqual); finalqual.append(overlapqual); finalqual.append(backqual);
		}
		//cout<<head<<endl<<ca<<endl<<cb<<endl<<refseq_left_use<<endl<<refseq_right_use<<endl;
		//cout<<head<<endl<<sequse<<endl<<seqcompare<<endl<<ca<<endl<<cb<<endl<<refseq_left_use<<endl<<refseq_right_use<<endl<<finalread<<endl<<endl;

	}}
	
	else{
		int flag2,pos2; string chr2,cigar2, tmp, tmp2;
		if(chr1a==chr1b and chr1a!=chr2a and chr1a==chr2b) // treat as 3-2 case, and fix read 1 using part of read 2	
			{flag2=flag2b; pos2=pos2b; chr2=chr2b; cigar2=cigar2b;}
		else if(chr1a==chr1b and chr1a!=chr2b and chr1a==chr2a) // treat as 3-2 case, and fix read 1 using part of read 2
			{flag2=flag2a; pos2=pos2a; chr2=chr2a; cigar2=cigar2a;}
		else if(chr2a==chr2b and chr2a!=chr1a and chr2a==chr1b) // treat as 2-3. requires switching read labels, for convenience.
			{
			flag2=flag1b; pos2=pos1b; chr2=chr1b; cigar2=cigar1b;
			flag1a=flag2a; flag1b=flag2b; pos1a=pos2a; pos1b=pos2b;
			chr1a=chr2a; chr1b=chr2b; string cigar1a=cigar2a; cigar1b=cigar2b;
			tmp=seq1; tmp2=qual1; seq1=seq2; qual1=qual2; seq2=tmp; qual2=tmp2;
			}
		else if(chr2a==chr2b and chr2a!=chr1b and chr2a==chr1a) // treat as 2-3. requires switching read labels, for convenience.
                        {
			flag2=flag1a; pos2=pos1a; chr2=chr1a; cigar2=cigar1a;
                        flag1a=flag2a; flag1b=flag2b; pos1a=pos2a; pos1b=pos2b;
                        chr1a=chr2a; chr1b=chr2b; cigar1a=cigar2a; cigar1b=cigar2b;
                        tmp=seq1; tmp2=qual1; seq1=seq2; qual1=qual2; seq2=tmp; qual2=tmp2;
			//cout<<seq2<<endl;
                        }
		else {output<<"@"<<head<<"\t"<<errorcase<<endl<<finalread<<endl<<"+"<<endl<<finalqual<<endl; continue;}
                vector<char> ciga_let, cigb_let;
                vector<int> ciga_num, cigb_num;		
		int matchpos, refpos; string refseq=""; // starting read position of interest, reference position relative to overlap, and reference sequence
		
		if(abs(pos1a-pos2) < abs(pos1b-pos2)) {mpos=pos1a; mflag=flag1a; mcigar=cigar1a; mpos=pos1a;}
                else {mpos=pos1b; mflag=flag1b; mcigar=cigar1b; mpos=pos1b;}
		
		if (mflag!=flag2 and abs(pos2-mpos)<=200){
		errorcase="2/2fixby1";
		parse_cigar(mcigar,ciga_num,ciga_let);  parse_cigar(cigar2,cigb_num,cigb_let);  // parse cigar strings
		get_matchseq(seq1,qual1,chr1a,mpos,ciga_num, ciga_let,fastapath,refseq,lclip); // get reference sequence
		runsw(1, 4, seq1, seq2, qual1, qual2, 250, ca, cb, qa, qb, starta, startb,olen); // get overlap
		
		// check if we can get additional nt from the ends of the string
		 string sloverhang=""; string  qloverhang=""; string sroverhang=""; string qroverhang="";
		if(startb>starta){
			for(int i=0; i<(startb-starta);i++)
			{
				sloverhang.push_back(seq2[i]); qloverhang.push_back(qual2[i]);	
			}
		}
		if((seq2.length()-startb)>(seq1.length()-starta))
		{
                        for(int i=seq1.length();i<seq2.length();i++)
                        {
                                sroverhang.push_back(seq2[i]); qroverhang.push_back(qual2[i]);
                        }
		}
		
		string frontpart, frontqual, backpart,backqual; // nonoverlapping parts
		for(int i=0;i<=starta-olen;i++){
			frontpart.push_back(seq1[i]);
			frontqual.push_back(qual1[i]);
		}
		for(int i=starta+1;i<seq1.length();i++){
			backpart.push_back(seq1[i]);
			backqual.push_back(qual1[i]);
		}

		// overlap
		int j=0; // keep track of reference nts
		int k=startb-olen+1; // overlapping indices for sequence 2
		string overlappart, overlapqual;
		for(int i=starta-olen+1; i<=starta; i++){
			if(seq1[i]==seq2[k] ) {overlappart.push_back(seq1[i]); overlapqual.push_back(qual1[i]);}
			else
			{
				if(i<lclip or i>=(lclip+refseq.length()))
				{ // where there is no reference, go by quality
                          		if(int(qual1[i])-int(qual2[k])>10){overlappart.push_back(seq1[i]); overlapqual.push_back(qual1[i]);}
                          		else if (int(qual2[k])-int(qual1[i])>10) {overlappart.push_back(seq2[k]); overlapqual.push_back(qual2[k]);}
                          		else {overlappart.push_back(seq1[i]); overlapqual.push_back(qual1[i]);}				
				}
				else { // use reference where possible
					j=i-lclip; // place in reference
					if (seq1[i]==refseq[j]) {overlappart.push_back(seq1[i]); overlapqual.push_back(qual1[i]);} // use matching reference
					else if (seq2[k]==refseq[j]) {overlappart.push_back(seq2[k]); overlapqual.push_back(qual2[k]);} // use matching reference
					else {if(seq1[i]=='N'){overlappart.push_back(seq2[k]); overlapqual.push_back(qual2[k]);} else {overlappart.push_back(seq1[i]); overlapqual.push_back(qual1[i]);}} // where neither matches reference use sequence 1, unless sequence 1 is N
				}
			}
			k++;
		}
		finalread.append(sloverhang); finalread.append(frontpart); finalread.append(overlappart); finalread.append(backpart); finalread.append(sroverhang);
		finalqual.append(qloverhang); finalqual.append(frontqual); finalqual.append(overlapqual); finalqual.append(backqual); finalqual.append(qroverhang);
		//cout<<head<<endl<<seq1<<endl<<seq2<<endl<<refseq<<endl<<ca<<endl<<cb<<endl<<finalread<<endl<<endl;	
		}
	}}

//if(errorcase.length()==0)
//	cout<<head<<endl;
//cout<<head<<endl<<errorcase<<endl;
if(finalread.length()>0)
	output<<"@"<<head<<"\t"<<errorcase<<endl<<finalread<<endl<<"+"<<endl<<finalqual<<endl;


}
//cout<<count<<endl;
readset1.close();
readset2.close();
output.close();

}

/**********************************************/
// other funcs

// tokenize the string
void tokenize(string read, vector<string> &read_vec, string delim){

    bool keepEmpty = false;
    
    size_t prev = 0;
    size_t next = 0;	

    while ((next = read.find_first_of(delim, prev)) != string::npos)
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

/*
string reverse(string seq, int seq_len){
        string revseq;
        char c;
       for (int i=seq_len-1; i>=0; i--){
           c=seq[i];
        revseq.push_back(c);
        }
        return revseq;
}
*/

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

void parse_cigar(string cigar, vector<int> &cigar_num, vector<char> &cigar_let ){


char *c = new char[cigar.length() + 1];
strcpy(c, cigar.c_str());
string val;

	for(int i=0;i<cigar.length();i++){
	   if(isdigit(c[i])) {val.push_back(c[i]);}
	   else {cigar_let.push_back(c[i]); cigar_num.push_back(atoi(val.c_str())); val="";}
	}

}

void get_matchseq(string seq,string qual, string chr, int pos, vector<int> cigar_num, vector<char> cigar_let,string fastapath, string &refseq, int &lclip){
// parse CIGAR string and grab reference sequence


// get clipping size and match length for initial match
	int i;
	for(i=0; i<cigar_let.size(); i++)
		if(cigar_let[i]=='M')
		{
			break;
		}
		else
			lclip=lclip+cigar_num[i];
			 //cout<<cigar_let[i]<<"\t"<<lclip<<"\t";

	// get reference fasta for matches
	for(int j=i;j<cigar_let.size(); j++){

		if(cigar_let[j]=='M')
		{
			refseq.append(searchfasta(chr, pos, pos+cigar_num[j]-1, fastapath)); // get reference (see next function)
			pos=pos+cigar_num[j];
		}
		else if (cigar_let[j]=='I')
		{
			string tmp;
			for(int k=0;k<cigar_num[j];k++)
				tmp.push_back('N');
			refseq.append(tmp);
		}
		else if (cigar_let[j]=='D')
			pos=pos+cigar_num[j];
	}	

	//cout<<i<<cigar_let[i]<<cigar_num[i]<<endl;

}


string searchfasta(string chr, int pos1,int pos2, string fastapath){
	// search the reference fasta file to grab the sequence for chr1:pos1-pos2 using samtools
	
	char seq[1000]; string rseq; ostringstream p1,p2;
	p1<<pos1; p2<<pos2;

	// write samtools command and convert to char
	string line="samtools faidx "+fastapath+" "+"\""+chr+"\"" +":"+p1.str()+"-"+p2.str().c_str()+"|sed 1d | tr -d \"\n\"";
	char *command = new char[line.length() + 1];
	strcpy(command, line.c_str());
	
	// run command
	FILE *fp;
	fp=popen(command,"r");
	fgets(seq, 1000, fp);
	pclose(fp);
	rseq=string(seq);
	
return rseq;
}
