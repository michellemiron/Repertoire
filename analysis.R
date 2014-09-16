## analyze TCR repertoire
# hypotheses: 
#
##1. After a certain period (6 months or later), the expanded T cell clones from donor stimulation are eliminated or reduced
##  alternatively, top clones in pre-tx stimulated condition are more likely to disappear than top clones in pre-tx unstimulated condition
##  Key caveat: sampling issue due to limited number of cells and sequencing reads. 
## solution: contorl for sampling issue by define "presence" using a threshold where the power of detection is large.

##2. The abundant clones after 6 months or later are not expanded after stimulation


args<-commandArgs(TRUE)

# cat(length(args), "\n")

file = args[1]
if (length(args)< 3 ) {
	freq1 = 1e-4  # freq1 = 1e-4,  "top"
	freq2 = 1e-5  # freq2 = 1e-5, "present"
} else {
	freq1 = as.numeric(args[2])
	freq2 = as.numeric(args[3])
}

tcr = read.table(file, header=T, sep="\t")

num = ncol(tcr)

## example: 
# nucleotide	total	ITN2_PreCD8_031313	ITNp2PreStim_CD8	ITN2_6moCD8_031313	ITN2_1yCD8_031313	ITN2_2yCD8_031313
# GGAGCTGGGGGACTCGGCCCTTTATCTTTGCGCCAGCAGCTTGGGGCCGGGTGGGGAGCA	1329732	0.024543792	61.76351143	0	0	0


# first summary stats
# number of clones observed
nclones = c()
for (i in 3:num) {
	nclones = c(nclones, length(tcr[tcr[,i]>0, 1]))
	
}
#cat("number of clones\n")
#cat(nclones, "\n")

#min observed freq
minfreq = c()
for (i in 3:num) {
	minfreq = c(minfreq, min(tcr[tcr[,i]>0, i])/100) 
}
#cat("min observed freq\n")
#cat(minfreq,"\n")	
	
# max freq

maxfreq = c()
for (i in 3:num) {
	maxfreq = c(maxfreq, signif(max(tcr[,i])/100, 3) )
}
#cat("max observed freq\n")
#cat(maxfreq, "\n")	


stats = t(rbind(nclones, minfreq, maxfreq))

print(stats, sep="\t")

## compare all clones between stimulated and unstimulated conditions: which is more likely to be present in post?

cat("test presence\n")
for (i in 5:num) {
# pre-tx unstimulated	
	x = c(nrow(tcr[tcr[,3] >= freq1 * 100 & tcr[,i] >= freq2 * 100 , ]), nrow(tcr[tcr[,3] >= freq1 * 100 & tcr[,i] < freq2 * 100 , ]))
	
# pre-tx stimulated
	y = c(nrow(tcr[tcr[,4] >= freq1 * 100 & tcr[,i] >= freq2 * 100 , ]), nrow(tcr[tcr[,4] >= freq1 * 100 & tcr[,i] < freq2 * 100 , ]))
	test = fisher.test(cbind(x,y))
#	print(x)
	cat(x,y, signif(test$p.value,2),  round(test$estimate,2),"\n")
}


## test 2: look at top expanded clones in stimulated condition, test if they are reduced in post
cat("Test if top expanded clones are reduced even than pre-tx unsti\n")

topExp = tcr[tcr[,4] >= freq1 * 100 & tcr[, 3] < tcr[,4] / 10 ,]

pretx = c(nrow(topExp[topExp[,3] >= freq2 * 100, ]), nrow(topExp[topExp[,3] < freq2 * 100, ]))

cat(pretx, "\n")
for (i in 5:num) {
	y = c(nrow(topExp[topExp[,i] >= freq2 * 100 , ]), nrow(topExp[topExp[,i] < freq2 * 100, ]))
	test = fisher.test(cbind(pretx, y))
	cat(y,signif(test$p.value,2),  round(test$estimate,2), "\n" )
}