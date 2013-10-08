r=$1
num=$2

qsub -N TCR8_$1 -cwd -l mem=10G,time=5:: -o log/ -e log/ downsample.sh TCR8_S1_L001_R1_001.trimed.bam $r $num
qsub -N TCR9_$1 -cwd -l mem=10G,time=5:: -o log/ -e log/ downsample.sh TCR9_S1_L001_R1_001.trimed.bam $r $num
qsub -N TCR23_$1 -cwd -l mem=10G,time=5:: -o log/ -e log/ downsample.sh TCR23_both.trimed.bam $r  $num
qsub -N TCR22_$1 -cwd -l mem=10G,time=5:: -o log/ -e log/ downsample.sh TCR22_both.trimed.bam $r  $num
qsub -N TCR7_$1 -cwd -l mem=10G,time=5:: -o log/ -e log/ downsample.sh TCR7_S1_L001_R1_001.trimed.bam $r $num
qsub -N TCR6_$1 -cwd -l mem=10G,time=5:: -o log/ -e log/ downsample.sh TCR6_S1_L001_R1_001.trimed.bam $r $num

