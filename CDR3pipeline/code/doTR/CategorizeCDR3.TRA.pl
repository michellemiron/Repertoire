###categorizing CDR3 by V J cassette
#!/usr/bin/perl
#use Bio::DB::Bam;
use strict;
use warnings;
#use Gtk2::Unique

open(Fin, './files/doTR/ReadCDR/TRAV.ref.IGMT.motif2.txt');
my @vref;
<Fin>; #remove the header line
while (<Fin>){
	chomp;
	#print "$_\n";
	push @vref, [split('\t', $_)]; # put trav name and 3'-end reference seq in two-dimensional array @ref3end
}
close(Fin);
open (FileNameJref,"./files/doTR/ReadCDR/TRAJ.ref.IGMT.motif3.txt");
my @jref;
<FileNameJref>; #remove the head line  
while (<FileNameJref>){
	chomp;
	push @jref, [split('\t', $_)];
	
}
close(FileNameJref);

my $datadir=$ARGV[0];
my $inFile = $ARGV[1];  #input file is filtered sam file

open (FileName, "$datadir/$inFile");

my %cdr3count;

#open (Fout1, ">$inFile.CDR3Length.txt");
while (<FileName>){
	chomp;
	#print "$_\n";
	my @cdrProt = split('\t', $_);
	if ($cdrProt[11]){
		#if($cdrProt[1] eq $cdrProt[2]){
			#print "$cdrProt[0]\t$cdrProt[1]\n";
		#}
		#if($vref[$cdrProt[3]][0] eq $jref[$cdrProt[4]][0]){
		#	print "$cdrProt[0]\t$$#vref[$cdrProt[3]]\n";
		#}
		#if ($cdrProt[11] ne "Out_of_frame"){
		$cdr3count{"byAllIndel.$cdrProt[3]_$cdrProt[4]_$cdrProt[11]_$cdrProt[17]_$cdrProt[18]_$cdrProt[19]"}++;
		$cdr3count{"byAll.$cdrProt[3]_$cdrProt[4]_$cdrProt[11]"}++;
		$cdr3count{"byCDR3.$cdrProt[11]"}++;
		$cdr3count{"byVJcomb.$cdrProt[3]_$cdrProt[4]"}++;
		$cdr3count{"byVsum.$cdrProt[3]"}++;
		$cdr3count{"byJsum.$cdrProt[4]"}++;
		#}
	}
}


close(FileName);
#close(Fout1);
#open(Fout3, ">$inFile.byVJ.txt");
open(Fout3, ">$datadir/$inFile.byEverythingNew.txt");
######print "\nGRADES IN ASCENDING NUMERIC ORDER:\n";
#foreach $key (sort hashValueAscendingNum (keys(%grades))) {
#   print "\t$grades{$key} \t\t $key\n";
#}

####print "\nGRADES IN DESCENDING NUMERIC ORDER:\n";
foreach my $key (sort hashValueDescendingNum (keys(%cdr3count))) {
   #print Fout3 "$cdr3count{$key}\t$key\n";
   my @Keycdr = split('\.',$key);
   my @output = split('_',$Keycdr[1]);
  
   print Fout3 "$cdr3count{$key}\t$Keycdr[0]\t";
	for (my $i=0; $i<=$#output;$i++){
   	print Fout3 "$output[$i]\t";
   }
   print Fout3 "\n";
}



sub hashValueDescendingNum {
  $cdr3count{$b} <=> $cdr3count{$a};
}
sub hashValueAscendingNum {
   $cdr3count{$a} <=> $cdr3count{$b};
}

close(Fout3);
