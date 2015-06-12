###########to find -vvvvv-jjjjj----pattern reads only in Yaping.FindVDJ6_4.AssignProb_CDR3_1.pl
###########To find both -vvvvvv---------- and -vvvvv-jjjjj------ pattern reads in this code: Yaping.FindVDJ6_4.AssignProb_CDR3_2.pl
###########The code originally named as "Yaping.FindVDJ6_4.AssignProb_CDR3_2.pl" for Mouse TRA. Adapted for human TRB
##########The code is modified from Yaping.ReadSam.CDR3DNAseqTRB_1.pl  In order to develop the bash pipeline, the input method is changed.
##########update @Vanno, @Vanno_ex, @Janno, @Janno_ex
##########v4: update to find the indels and mismatchs of TRBV and TRBJ
##########v5: Fix the problem that productive CDR3 sequences lose 
##########v6: add V gene mapping quality score
#!/usr/bin/perl
#use Bio::DB::Bam;
use strict;
use warnings;
#use List::MoreUtils qw(uniq);
#use Gtk2::Unique


### Auto_120305: 2792192 reads; Auto_..87: 1584373 reads; Auto_..90: 1504752 reads; 
### Auto_120305: Auto_..87 : Auto..90 = 1.86:1.05:1


my @Readsname;
my $NumFile=0;


open(FinIgHpos,'./files/doTR/ReadSam/TRB.v37p8.YapingCoordinates2.txt');
my $temps = <FinIgHpos>; # remove the header line
my @IgHanno;
while (<FinIgHpos>) {
	chomp;
	#print "$_\n";
	push @IgHanno, [split('\t', $_)]; # put annotation in two-dimensional array @IgHanno
	
}
close(FinIgHpos);



#my $outdir='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcell/dat/processedbams/human/TCRB';

my $coreName;
my $datadir=$ARGV[0];
my $inFile = $ARGV[1];  #input file is filtered sam file

open (Fout, ">$datadir/$inFile.CDR3.VJ.seq");
open (FoutV, ">$datadir/$inFile.CDR3.Vonly.seq");
open (FoutJ, ">$datadir/$inFile.CDR3.Jonly.seq");
open (FoutJ4cm, ">$datadir/$inFile.Jonly.4cm.fasta");
open (mF, "$datadir/$inFile");

my $indxP=0;
my @cdrV_ex;
my @cdrJ_ex;
my @cdrV;
my @cdrJ;
my $ReadsID_ex;
my @VDJcount;
my @VPos;
my @JPos;
my @values_ex;
my @VJanno_ex;
my @Vanno_ex;
my @Vmscore_ex;
my @Vmscore;
my @Janno_ex;
my @VJanno;
my @Vanno;
my @Janno;

#my $refSeqV;
#my $readSeqV;
#my $refSeqJ;
#my $refSeqJ;
#my $refSeqV_ex;
#my $readSeqV_ex;
#my $refSeqJ_ex;
#my $readSeqJ_ex;
my $refSeq;
my $readSeq;
my $refSeq_ex;
my $readSeq_ex;
my @VdJdI = (-10, -10, -10, -10, -5, -10, -5,  -10, -5, -10, -5);
my @VdJdI_ex;
##this array stores # of V deletion [0], # of J deletion [1], # of insertion [2], V frame del [3], V frame del pos [4], 
##V frame insertion [5], V frame insertion pos[6],
## J frame del [7], J frame del pos [8], J frame insertion [9], J frame insertion pos [10].
my @VmJm = ('f','f','f','f');	# V frame mismatch [0], V frame mismatch pos [1],  J frame mismatch [2], J frame mismatch pos [3]						
my @VmJm_ex;

while (<mF>){
	
#print "$_\n";
	chomp;
	my @values = split('\t', $_);
#print "$i\t$values[13]\n";
	my @num; #for reference
	my @cigar; #for reference
	my @startT;  #start match position in reference
	my @endT;	#end match position in reference
	my @startTM;  #start match position in reference for 'M' cigar only
	my @endTM; #end match position in reference for 'M' cigar only
	my $indxJ=0;

	$startT[0] = $values[3];
	$endT[0] = $startT[0];
	
	my @startR;#start match postion in short reads 
	my @endR; #end match postion in short reads 
	my @startRM; #start match postion in short reads for 'M' cigar only
	my @endRM; #end match postion in short reads for 'M' cigar only
	my @numR; #for short reads
	my @cigarR; #for short reads
	
	my @vindelmm;
	my @jindelmm;

	if ($#values > 3 && $values[9] !~ /N/) { # remove SQ lines & sequences with "N"
		$startT[0] = $values[3];
		$endT[0] = $startT[0];	
		while ( $values[5]=~m/(\d+)([MDN])/g){ ##to find out the mapping position in reference genome --chromosome
		#while ($values[5]=~m/(\d+)([MIDNSHP=X])/g){  #for print cigar pattern
			push (@num, $1);  #number of match, deletion, insertion, ...
			push (@cigar, $2); ####MIDNSHP=X
			
			$endT[$indxJ] = $startT[$indxJ]+$1-1;
			$indxJ++;
			$startT[$indxJ]=$endT[$indxJ-1]+1;
			
			if ($cigar[$indxJ-1] eq 'M'){
				push(@endTM, $endT[$indxJ-1]);								
				push(@startTM, $startT[$indxJ-1]);	
			}
			
		}
		
		$indxJ = 0;
		$startR[0] = 1;
		$endR[0] = 0;
		while ($values[5]=~m/(\d+)([MIS])/g){
		#while ($values[5]=~m/(\d+)([MIDNSHP=X])/g){  #for print cigar pattern
			push (@numR, $1);  #number of match, deletion, insertion, ...
			push (@cigarR, $2); ####MIDNSHP=X
			
			$endR[$indxJ] = $startR[$indxJ]+$1-1;
			$indxJ++;
			$startR[$indxJ]=$endR[$indxJ-1]+1;
			if ($cigarR[$indxJ-1] eq 'M'){
				push(@endRM, $endR[$indxJ-1]);								
				push(@startRM, $startR[$indxJ-1]);	
			}
			
		}
		#print "$values[0]\t$values[5]\t";
		#print join ("*",@cigarR);
		#print "--";
		#print join("*",@cigar);
		#print "\n";
		#print join (" ",@startRM);
		#print "--";
		#print join (" ", @startTM);
		#print "\n";
		#print join ("_",@endRM);
		#print "--";
		#print join ("_", @endTM);
		#print "\n";
				
		for (my $i=0;$i<=$#startRM;$i++) {
				for (my $j=0; $j<$#IgHanno; $j++){
					if (($startTM[$i]>=$IgHanno[$j][2]&&$startTM[$i]<=$IgHanno[$j][3])||($endTM[$i]>=$IgHanno[$j][2]&&$endTM[$i]<=$IgHanno[$j][3])){
						
						if (substr($IgHanno[$j][1],0,4) eq "TRBV"){  #$IgHanno[][7] is the DNA sequence for the gene
							push(@VPos,$values[1]);		
							push(@VPos, $startRM[$i]);
							push(@VPos, $endRM[$i]);
							push (@Vanno, $IgHanno[$j][1]);
							push(@Vmscore, $values[4]);
							
							last;
						}	
						elsif (substr($IgHanno[$j][1],0,4) eq "TRBJ"){
							push(@JPos,$values[1]);
							push(@JPos, $startRM[$i]);
							push(@JPos, $endRM[$i]);
							push (@Janno, $IgHanno[$j][1]);
							last;
						}
						#print "Start\t$values[0]\t$IgHanno[$j][1]\t$startM[$i]\t"; 
					}		
				}  
		}
		############### to find indels and mismatches --START ###########################
		my $tmpSpos;
		my $tmpEpos;
		for (my $i=0;$i<=0;$i++) { ##just consider the first M cigar
				for (my $j=0; $j<$#IgHanno; $j++){
					#print "$startTM[$i]\n";
                                        #print "$endTM[$i]\n";
                                        #print "$IgHanno[$j][3]\n";
					#sleep(10);
					if (($startTM[$i]>=$IgHanno[$j][2]&&$startTM[$i]<=$IgHanno[$j][3])||($endTM[$i]>=$IgHanno[$j][2]&&$endTM[$i]<=$IgHanno[$j][3])){
						
						if (substr($IgHanno[$j][1],0,4) eq "TRBV"){  #$IgHanno[][7] is the DNA sequence for the gene
							$tmpSpos=length($IgHanno[$j][7])-($IgHanno[$j][3]-$values[3])-1;
							$tmpEpos=length($IgHanno[$j][7])-$tmpSpos+1;
							$refSeq = substr($IgHanno[$j][7],$tmpSpos,$tmpEpos);
							if ($tmpSpos<0){
								$refSeq = substr($IgHanno[$j][7],0,$tmpEpos);
							}
							else{
								$refSeq = substr($IgHanno[$j][7],$tmpSpos,$tmpEpos);
							}
							$readSeq = substr($values[9], $startRM[0]-1,$endRM[$#endRM]-$startRM[0]+1);
							my $numVd = length($refSeq)-length($readSeq);
							if ($numVd < 0){
							    $VdJdI[0] = 0; ##this array stores # of V deletion, # of J deletion, # of insertion, V cigar, J cigar.
							}
							else {
								$VdJdI[0] = $numVd; 
							}
							$indxJ = 0;
							my @tmpstartR;
							my @tmpendR;
							$tmpstartR[0] = 1;
							$tmpendR[0] = 0;
							##this array stores # of V deletion [0], # of J deletion [1], # of insertion [2], V frame del [3], V frame del pos [4], 
							##V frame insertion [5], V frame insertion pos[6], 
							## J frame del [7], J frame del pos [8], J frame insertion [9], J frame insertion pos [10].
							while ($values[5] =~ m/(\d+)([MID])/g){
								$tmpendR[$indxJ] = $tmpstartR[$indxJ]+$1-1;
								$indxJ++;
								#print "$1\t$2\n";
								$tmpstartR[$indxJ]=$tmpendR[$indxJ-1]+1;
								if ($2 eq 'D'){
									$VdJdI[3] = $1;
									$VdJdI[4] = $tmpstartR[$indxJ-1];
									if (length($refSeq)>=($VdJdI[4]-1)){
									for (my $k=0; $k<$1; $k++){
										substr($readSeq, $VdJdI[4]-1, 0) = '*';								
									}
									}
								}
								if ($2 eq 'I'){  #75S9M1I45M85S
									$VdJdI[5] = $1;
									$VdJdI[6] = $tmpstartR[$indxJ-1];
									#print "$VdJdI[5]\t$VdJdI[6]***$refSeq\n";
									if (length($refSeq)>=($VdJdI[6]-1)){
									for (my $k=0; $k<$1; $k++){
										substr($refSeq, $VdJdI[6]-1, 0) = '*';	 #78S59M2I8M78S 	73S61M1I8M77S						
									}
									}								
								}
								
							}
							my $tmpMMpos = "f";
							my $tmpMMch = "f";
							for (my $k=0; $k<length($refSeq); $k++){
								if ($k < length($readSeq)){
								my $tmpGch = substr($refSeq, $k, 1);
								my $tmpRch = substr($readSeq, $k, 1);
								#$tmpRch = uc $tmpGch;
								if (uc $tmpGch ne $tmpRch  && $tmpGch ne '*' && $tmpRch ne '*'){
									$tmpMMch.="$tmpGch:$tmpRch;";
									$tmpMMpos.="$k;";
									$tmpMMch =~ s/f//;
									$tmpMMpos =~ s/f//;
									
								}
								}
							}
							$VmJm[0] = $tmpMMch;
							$VmJm[1] = $tmpMMpos;
							last;
						}	
						elsif (substr($IgHanno[$j][1],0,4) eq "TRBJ"){
							$tmpSpos=length($IgHanno[$j][7])-($IgHanno[$j][3]-$values[3])-1;
							$tmpEpos=length($IgHanno[$j][7])-$tmpSpos+1;
							$refSeq = substr($IgHanno[$j][7],$tmpSpos,$tmpEpos);
							if ($tmpSpos<0){
								$refSeq = substr($IgHanno[$j][7],0,$tmpEpos);
							}
							else{
								$refSeq = substr($IgHanno[$j][7],$tmpSpos,$tmpEpos);
							}
							$readSeq = substr($values[9], $startRM[0]-1,$endRM[$#endRM]-$startRM[0]+1);
							
							if ($tmpSpos<0){
								$VdJdI[1] = 0;
							}
							else {
								$VdJdI[1] = $tmpSpos;
							}
							$indxJ = 0;
							my @tmpstartR;
							my @tmpendR;
							$tmpstartR[0] = 1;
							$tmpendR[0] = 0;
							##this array stores # of V deletion [0], # of J deletion [1], # of insertion [2], V frame del [3], V frame del pos [4], 
							##V frame insertion [5], V frame insertion pos[6], 
							## J frame del [7], J frame del pos [8], J frame insertion [9], J frame insertion pos [10].
							while ($values[5] =~ m/(\d+)([MID])/g){
								$tmpendR[$indxJ] = $tmpstartR[$indxJ]+$1-1;
								$indxJ++;
								$tmpstartR[$indxJ]=$tmpendR[$indxJ-1]+1;
								if ($2 eq 'D'){
									$VdJdI[7] = $1;
									$VdJdI[8] = $tmpstartR[$indxJ-1];
									if (length($refSeq)>=($VdJdI[8]-1)){
									for (my $k=0; $k<$1; $k++){
										substr($readSeq, $VdJdI[8]-1, 0) = '*';								
									}
									}
								}
								if ($2 eq 'I'){
									$VdJdI[9] = $1;
									$VdJdI[10] = $tmpstartR[$indxJ-1];
									if (length($refSeq)>=($VdJdI[10]-1)){
									for (my $k=0; $k<$1; $k++){
										substr($refSeq, $VdJdI[10]-1, 0) = '*';								
									}
									}								
								}
								
							}
							my $tmpMMpos = "f";
							my $tmpMMch = "f";
							
							
							for (my $k = 0; $k<length($refSeq); $k++){
								if ($k < length($readSeq)){							
								my $tmpGch = substr($refSeq, $k, 1);
								my $tmpRch;
								if ($tmpSpos < 0){
									my $tmpk = $k+abs($tmpSpos);
									$tmpRch = substr($readSeq, $tmpk, 1);
								}
								else {
									$tmpRch = substr($readSeq, $k, 1);
								}
								#$tmpRch = uc $tmpGch;
								if (uc $tmpRch ne $tmpRch && $tmpGch ne '*' && $tmpRch ne '*'){
									$tmpMMch.="$tmpGch:$tmpRch;";
									$tmpMMpos.="$k;";
									$tmpMMch =~ s/f//;
									$tmpMMpos =~ s/f//;
									
								}
								}
							}
							$VmJm[2] = $tmpMMch;
							$VmJm[3] = $tmpMMpos;
							last;
						}
						#print "Start\t$values[0]\t$IgHanno[$j][1]\t$startM[$i]\t"; 
					}		
				}  
				
		}
		#print join ("vv",@VPos);
		#print "\t";
		#print join ("jj",@JPos);
		#print "\n";
		#print "$values[0]\t$tmpSpos\t$values[5]\t";
		#print join ("_", @VdJdI);
		#print "\t";
		#print join ("_", @VmJm);
		
		#print "\n$refSeq\n$readSeq\n";
		############### to find indels and mismatches--END ############################
		
		if ($indxP==0){
			 @cdrV = @VPos;
			 @cdrJ = @JPos;
			 $ReadsID_ex = $values[0];
			 @values_ex = @values;
			 @VJanno_ex = @VJanno;
			 @Janno_ex = @Janno;
			 @Vanno_ex = @Vanno;
			 @VdJdI_ex = @VdJdI;
			 @VmJm_ex = @VmJm;
			 @Vmscore_ex = @Vmscore;
		}
		elsif ($values[0] eq $ReadsID_ex) {
			@cdrV = (@cdrV, @VPos);
			@cdrJ = (@cdrJ, @JPos);
			@VJanno_ex = (@VJanno_ex, @VJanno);
			@Vanno_ex = (@Vanno_ex, @Vanno);
			@Janno_ex = (@Janno_ex, @Janno);
			$ReadsID_ex = $values[0];
			@values_ex = @values;
			@Vmscore_ex =(@Vmscore_ex, @Vmscore);
			for (my $k=0; $k<=$#VdJdI; $k++){
				if ($VdJdI_ex[$k]==-10 || $VdJdI_ex[$k]==-5){
					$VdJdI_ex[$k] = $VdJdI[$k]; 
				}
			}
			for (my $k=0; $k<=$#VmJm; $k++){
				if ($VmJm_ex[$k] eq 'f' ){
					$VmJm_ex[$k] = $VmJm[$k]; 
				}
			}
			#print join ("_", @VmJm);
			#print "\n";
		}
		elsif ($values[0] ne $ReadsID_ex){

			#print Fout1 "$ReadsID_ex\t";
			#print Fout1 "@IgHuniq_ex";	
			#print Fout1 "\n";
			my $cdrStart;
			my $cdrEnd;
			my $cdrDna;
			if ($#cdrV>=2 && $#cdrJ>=2){
				if ($cdrV[0]==16 && $cdrJ[0]==16){
					$cdrStart = $cdrV[1];
					$cdrEnd = $cdrJ[2];
					$VdJdI_ex[2] = $cdrJ[1]-$cdrV[$#cdrV]-1;
					for (my $i=3; $i<=$#cdrV; $i+=3){
						if ($cdrStart > $cdrV[$i+1]){
							$cdrStart = $cdrV[$i+1];
						}	
					}
					for (my $i=3; $i<=$#cdrJ; $i+=3){
						if ($cdrEnd < $cdrJ[$i+2]){
							$cdrEnd = $cdrJ[$i+2];
						}	
					}
					
					if ($values_ex[1]==16){
						$cdrDna = substr($values_ex[9],$cdrStart-1, $cdrEnd-$cdrStart+1);
						print Fout "$ReadsID_ex\t16\t";
						print Fout join (".", @Vanno_ex);
						print Fout ".";
						print Fout join (".", @Janno_ex);
						print Fout "\t$cdrDna\t";
						print Fout "$cdrStart\t$cdrV[$#cdrV]\t$cdrJ[1]\t$cdrEnd\t";						
						print Fout join ("\t", @VdJdI_ex);
						print Fout "\t";
						print Fout join ("\t", @VmJm_ex);
						print Fout "\t";
						print Fout join ("_", @Vmscore_ex);
						print Fout "\n";
						
					}
				}
				if ($cdrV[0]==0 && $cdrJ[0]==0){
					$cdrStart = $cdrV[1];
					$cdrEnd = $cdrJ[2];
					$VdJdI_ex[2] = $cdrJ[1]-$cdrV[$#cdrV]-1;					
					for (my $i=3; $i<=$#cdrV; $i+=3){
						if ($cdrStart < $cdrV[$i+1]){
							$cdrStart = $cdrV[$i+1];
						}	
					}
					for (my $i=3; $i<=$#cdrV; $i+=3){
						if ($cdrEnd > $cdrJ[$i+2]){
							$cdrEnd = $cdrJ[$i+2];
						}	
					}
					if ($values_ex[1]==0){
						$cdrDna = substr(reverseIUPAC($values_ex[9]),length($values_ex[9])-$cdrStart-1, $cdrStart-$cdrEnd+1);
						print Fout "$ReadsID_ex\t2\t";
						print Fout join (".", @Vanno_ex);
						print Fout ".";
						print Fout join (".", @Janno_ex);
						print Fout "\t$cdrDna\t";
						print Fout "$cdrStart\t$cdrV[$#cdrV]\t$cdrJ[1]\t$cdrEnd\t";
						print Fout join ("\t", @VdJdI_ex);
						print Fout "\t";
						print Fout join ("\t", @VmJm_ex);
						print Fout "\t";
						print Fout join ("_", @Vmscore_ex);
						print Fout "\n";
						
						
					}
					
				}
				
			}
			elsif ($#cdrV>=2 && $#cdrJ==-1){   #### V only
				if ($cdrV[0]==16 ){
					$cdrStart = $cdrV[1];
					$cdrEnd = $cdrV[$#cdrV]+70;
					for (my $i=3; $i<=$#cdrV; $i+=3){
						if ($cdrStart > $cdrV[$i+1]){
							$cdrStart = $cdrV[$i+1];
							$cdrEnd = $cdrV[$i+2]+70;
						}	
					}
					
					if ($values_ex[1]==16){  
						$cdrDna = substr($values_ex[9],$cdrStart-1, $cdrEnd-$cdrStart+1);
						print FoutV "$ReadsID_ex\t16\t";
						print FoutV join (".", @Vanno_ex);
						print FoutV "\t$cdrDna\t";
						print FoutV "$cdrStart\t$cdrEnd\t$values_ex[9]\t";
						print FoutV join ("_", @Vmscore_ex);
						print FoutV "\n";
					}
				}
			}
			elsif ($#cdrV==-1 && $#cdrJ>=2){  #### J only
				if ($cdrJ[0]==16 ){
		
					$cdrStart = $cdrJ[1]-60;
					$cdrEnd = $cdrJ[2];		
					for (my $i=3; $i<=$#cdrJ; $i+=3){
						if ($cdrEnd < $cdrJ[$i+2]){
							$cdrEnd = $cdrJ[$i+2];
							$cdrStart = $cdrJ[$i+1]-60;
						}	
					}
					if ($cdrStart < 0){
						$cdrStart = 1;
					}
					if ($values_ex[1]==16){
						$cdrDna = substr($values_ex[9],$cdrStart-1, $cdrEnd-$cdrStart+1);
						print FoutJ "$ReadsID_ex\t16\t";
						print FoutJ join (".", @Janno_ex);
						print FoutJ "\t$cdrDna\t";
						print FoutJ "$cdrStart\t$cdrEnd\t$values_ex[9]\n";
						print FoutJ4cm ">$ReadsID_ex\t16\t";
						print FoutJ4cm join (".", @Janno_ex);
						print FoutJ4cm "\n$cdrDna\n";
						
					}
				}
			}
			@cdrV = @VPos;
			@cdrJ = @JPos;
    		$ReadsID_ex = $values[0];
    		@values_ex = @values;
    		
    		$#VJanno_ex = -2;
    		$#Vanno_ex = -2;
    		$#Janno_ex = -2;
    		$#VdJdI_ex = -2;
    		$#Vmscore_ex = -2;
    		
    		@VJanno_ex = @VJanno;	
    		@Janno_ex = @Janno;
    		@Vanno_ex = @Vanno;
    		@VdJdI_ex = @VdJdI;
    		@VmJm_ex = @VmJm;
    		@Vmscore_ex = @Vmscore;
    	}
		$#VPos=-2;
		$#JPos=-2;
		$#VJanno = -2;
		$#Vanno = -2;
		$#Janno = -2;
		$#Vmscore = -2;
		@VdJdI = (-10, -10, -10, -10, -5, -10, -5, -10, -5, -10, -5);
		@VmJm = ('f','f','f','f');
		$indxP++;  #to count if it is the first line except SQ@ lines	

	} ###if ($#values>3)
	} ### while (..)
close(Fout);
close(FoutV);
close(FoutJ);
close(mF);	
close(FoutJ4cm);

sub reverseIUPAC {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/NABCDGHMNRSTUVWXYabcdghmnrstuvwxy/NTVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
}
