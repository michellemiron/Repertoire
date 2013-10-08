infile=$1
entropydir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/human/entropy  # $2
VJdir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/human/VJcounts  # $3
#entropydir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/human/matchentropy
#VJdir=/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/human/matchVJcounts



#CDR3="{'`cut -f5 $infile |perl -pe 'chomp if eof'| sed 1d| tr \"\n\" \" \" | sed \"s/ /\' \'/g\"`'}"
#VJ="{'`cut -f1 $infile |perl -pe 'chomp if eof'| sed 1d| tr \"\n\" \" \" | sed \"s/ /\' \'/g\"`'}"
#V="{'`cut -f3 $infile |perl -pe 'chomp if eof'| sed 1d| tr \"\n\" \" \" | sed \"s/ /\' \'/g\"`'}"
#J="{'`cut -f4 $infile |perl -pe 'chomp if eof'| sed 1d| tr \"\n\" \" \" | sed \"s/ /\' \'/g\"`'}"
#counts="[`cut -f2 $infile | sed 1d| tr \"\n\" \" \" `]"

mkdir -p $entropydir
mkdir -p $VJdir
~/qsubmat -l mem=2G,time=:30: -e log -o log "get_entropy2(0,'$infile','$entropydir','$VJdir')"
~/qsubmat -l mem=2G,time=:30: -e log -o log "get_entropy2(1,'$infile','$entropydir','$VJdir')"
~/qsubmat -l mem=2G,time=:30: -e log -o log "get_entropy2(2,'$infile','$entropydir','$VJdir')"
~/qsubmat -l mem=2G,time=:30: -e log -o log "get_entropy2(3,'$infile','$entropydir','$VJdir')"
~/qsubmat -l mem=2G,time=:30: -e log -o log "get_entropy2(0,'$infile','$entropydir','$VJdir',1000)"
~/qsubmat -l mem=2G,time=:30: -e log -o log "get_entropy2(1,'$infile','$entropydir','$VJdir',1000)"
~/qsubmat -l mem=2G,time=:30: -e log -o log "get_entropy2(2,'$infile','$entropydir','$VJdir',1000)"
~/qsubmat -l mem=2G,time=:30: -e log -o log "get_entropy2(3,'$infile','$entropydir','$VJdir',1000)"

