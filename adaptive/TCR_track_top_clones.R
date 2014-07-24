## analyze TCR repertoire
# hypotheses: 
#
##1. After a certain period (6 months or later), the expanded T cell clones from donor stimulation are eliminated or reduced
##  alternatively, top clones in pre-tx stimulated condition are more likely to disappear than top clones in pre-tx unstimulated condition
##  Key caveat: sampling issue due to limited number of cells and sequencing reads. 
## solution: contorl for sampling issue by define "presence" using a threshold where the power of detection is large.

##2. The abundant clones after 6 months or later are not expanded after stimulation


compare <- function(tcr, freq1, freq2, fold=10) {
    if (ncol(tcr) < 2) {
        return(NA)
    }
                                        #	x = tcr[which(rowMeans(tcr) > 0.001),]
    x = tcr	
    x = normalize(x)
                                        #	print(dim(x))
    nsample = ncol(x)
    
#    stimCol <- grep("tim", colnames(tcr))
    stimCol <- grep("pretx_stim", colnames(tcr))
    
     preCol <- grep("pretx_unstim", colnames(tcr))[1]

 #   postCols <- grep("pretx", colnames(tcr), invert=T)
    # postCols <- grep("tim|tx", colnames(tcr), invert=T)
    
                                        # top clones from in vitro stimulated sample
    # topSti = tcr[tcr[,stimCol] > tcr[, preCol] * fold & tcr[,stimCol] >=  freq1, ]
### pre-stim defined donor reactive clones 
    topStiRows = (tcr[,stimCol[1]] > tcr[, preCol] * fold & tcr[,stimCol[1]] >=  freq1 )
    if (length(stimCol) > 1) {
        for (i in 2:length(stimCol)) {
            topStiRows = (topStiRows | (tcr[,stimCol[i]] > tcr[, preCol] * fold & tcr[,stimCol[i]] >=  freq1 ))
        }
    }
    
                                       
                         # topSti = tcr[tcr[,stimCol] > tcr[, preCol] * fold & tcr[,stimCol] >= freq1, ]
    topSti = tcr[topStiRows,]
    allPre = tcr[tcr[,preCol] >= freq2 , ]
    

### post-stim defined donor reactive clones
	postStimCol <- grep("posttx_antidonor", colnames(tcr))
	postUnstimCol <- grep("12mo_unstim", colnames(tcr))
	
	topPostStiRows <- (tcr[,postStimCol] > tcr[, postUnstimCol] * fold & tcr[,postStimCol] >=  freq1 )
	topPostSti = tcr[topPostStiRows, ]

	npreRx = nrow(topPostSti[topPostSti[,stimCol]> freq1, ])
	npostRx = nrow(topSti[topSti[, postStimCol]> freq1 , ])

	pdf("postReactive.pdf")
	maxFreq = max(max(tcr[,stimCol], tcr[,postStimCol]))
#	plot(c(topPostSti[1, stimCol], topPostSti[1, postStimCol]), xlab="Samples", ylab="Frequency", ylim=range(c(0,maxFreq)),  xaxt = "n", col='white')

	plot(c(log10(topPostSti[1, stimCol] + 0.0000001), log10(topPostSti[1, postStimCol] + 0.0000001)), xlab="Samples", ylab="Frequency", ylim=range(c(-7, 0)),  xaxt = "n", col='white')
	
	for (i in 1:nrow(topPostSti)) {
		points(c(log10(topPostSti[i, stimCol] + 0.0000001), log10(topPostSti[i, postStimCol] + 0.0000001)), pch=19, col=i)
		lines(c(log10(topPostSti[i, stimCol] + 0.0000001), log10(topPostSti[i, postStimCol] + 0.0000001)), lty="dotted", col=i)
		
	}

    axis(1, at=1:2, labels=c("pre-tx stim", "post-tx stim"))
	#text(1.5, 1.5, paste("npreRx: ", npreRx, "\nnpostRx: ", npostRx))
	dev.off()



	pdf("preReactive.pdf")
#	plot(c(topPostSti[1, stimCol], topPostSti[1, postStimCol]), xlab="Samples", ylab="Frequency", ylim=range(c(0,maxFreq)),  xaxt = "n", col='white')

	plot(c(log10(topSti[1, stimCol] + 0.0000001), log10(topSti[1, postStimCol] + 0.0000001)), xlab="Samples", ylab="Frequency", ylim=range(c(-7, 0)),  xaxt = "n", col='white')
	
	for (i in 1:nrow(topSti)) {
		points(c(log10(topSti[i, stimCol] + 0.0000001), log10(topSti[i, postStimCol] + 0.0000001)), pch=19, col=i)
		lines(c(log10(topSti[i, stimCol] + 0.0000001), log10(topSti[i, postStimCol] + 0.0000001)), lty="dotted", col=i)
		
	}

    axis(1, at=1:2, labels=c("pre-tx stim", "post-tx stim"))
	#text(1.5, 1.5, paste("npreRx: ", npreRx, "\nnpostRx: ", npostRx))
	dev.off()
}





normalize <- function(tcr){
    if (ncol(tcr) < 1) {
        return(NA)
    }
    for (f in 1:ncol(tcr)) {
		tcr[,f] = tcr[,f] / sum(tcr[,f]) 

	}
	return(tcr)	
}


# args<-commandArgs(TRUE)
library(getopt)

spec <- matrix(c(
        'input'     , 'i', 1, "character", "input csv or tsv file (required)",
        'type'     , 't', 2, "character", "csv or tsv (optional); default: tsv",
        "freq1" , "m", 2, "double", "min freq to be considered as top clones",
        "freq2", "d", 2, "double", "min freq to be considered as detectable", 
        "fold", "f", 2, "integer", "min fold change to be considered as expanded", 
        "vj", 'v', 2, "integer", "[1/0] do analysis on V/J (default 1 means yes)",
        'help'   , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)

opt = getopt(spec);

if (!is.null(opt$help) || is.null(opt$input)) {
    cat(paste(getopt(spec, usage=T),"\n"));
    q();
}

if (is.null(opt$vj)) {
	opt$vj = 0
}

# cat(length(args), "\n")

file = opt$input

freq1 = 1e-4  # freq1 = 1e-4,  "top"
freq2 = 1e-5  # freq2 = 1e-5, "present"
fold = 10  # fold change default threshold

itype = "t"
if (!is.null(opt$type)) {
	itype = opt$type
	
}

if(!is.null(opt$freq1)) {
	freq1 = opt$freq1
}

if(!is.null(opt$freq2)) {
	freq2 = opt$freq2
}

if(!is.null(opt$fold)) {
	fold = opt$fold
	
}

if( itype == "c")  {
	tcr = read.table(file, header=T, sep=",")
} else if (itype == "t") {
	tcr = read.table(file, header=T, sep="\t")			
} else {
	tcr = read.table(file, header=T)			
}

num = ncol(tcr)

## example: 
# nucleotide	total	ITN2_PreCD8_031313	ITNp2PreStim_CD8	ITN2_6moCD8_031313	ITN2_1yCD8_031313	ITN2_2yCD8_031313
# GGAGCTGGGGGACTCGGCCCTTTATCTTTGCGCCAGCAGCTTGGGGCCGGGTGGGGAGCA	1329732	0.024543792	61.76351143	0	0	0


cd4 <- tcr[, grep("CD4|cd4", colnames(tcr))]
cd8 <- tcr[, grep("CD8|cd8", colnames(tcr))]

cd4 <- normalize(cd4)
cd8 <- normalize(cd8)

# first summary stats
# cd4stats = repSum(cd4)
# cd8stats = repSum(cd8)

# print(cd4stats[order(rownames(cd4stats), decreasing=T), ])
# print(cd8stats[order(rownames(cd8stats), decreasing=T), ])

#if (!is.null(cd4) & dim(cd4)[2] > 1) {
compare(cd4, freq1, freq2, fold)
#}

#if(!is.null(cd8) & dim(cd8)[2] > 1) { 
#		 compare(cd8, freq1, freq2, fold)
# }

