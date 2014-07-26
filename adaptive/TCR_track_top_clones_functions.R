## analyze TCR repertoire
# hypotheses: 
#
##1. After a certain period (6 months or later), the expanded T cell clones from donor stimulation are eliminated or reduced
##  alternatively, top clones in pre-tx stimulated condition are more likely to disappear than top clones in pre-tx unstimulated condition
##  Key caveat: sampling issue due to limited number of cells and sequencing reads. 
## solution: contorl for sampling issue by define "presence" using a threshold where the power of detection is large.

##2. The abundant clones after 6 months or later are not expanded after stimulation


compare <- function(tcr, freq1, freq2, fold=5, prefix = "TopClone") {
    if (ncol(tcr) < 2) {
        return(NA)
    }

    nsample = ncol(tcr)

	# pre-tx stim sample
    stimCol <- grep("pretx_stim", colnames(tcr))
    
    # pre-tx unstim sample
    preCol <- grep("pretx_unstim", colnames(tcr))[1]

	# unstim samples
	unstimCols <- grep("tx_stim|tx_anti", colnames(tcr), invert=T)
	
	# post-tx stim sample in MLR
   	postStimCol <- grep("posttx_antidonor", colnames(tcr))
   	# post-tx unstim sample used in MLR
	postUnstimCol <- grep("12mo_unstim", colnames(tcr))


     # top clones from in vitro stimulated sample
    
### pre-tx stim defined donor reactive clones; row index
    topStiRows = (tcr[,stimCol[1]] > tcr[, preCol] * fold & tcr[,stimCol[1]] >=  freq1 )
    topSti = tcr[topStiRows,]
    

### post-tx stim defined donor reactive clones
	
	topPostStiRows <- (tcr[,postStimCol] > tcr[, postUnstimCol] * fold & tcr[,postStimCol] >=  freq1 )
	topPostSti = tcr[topPostStiRows, ]

### number of post-tx Stim donor reactive clones that are frequent in pre-tx stim sample
	npreRx = nrow(topPostSti[topPostSti[,stimCol]>= freq1, ])
	
### number of pre-tx Stim donor reactive clones that are frequent in post-tx stim sample
	npostRx = nrow(topSti[topSti[, postStimCol] >= freq1 , ])


### trend of reactive clones from post-tx stim
	pdf(paste(prefix, "postReactive.pdf", sep = "_"))
	maxFreq = max(max(tcr[,stimCol], tcr[,postStimCol]))
	totalFreqPost = sum(topPostSti[,postStimCol])
	totalFreqPre = sum(topPostSti[,stimCol])
## plot in original scale	
#	plot(c(topPostSti[1, stimCol], topPostSti[1, postStimCol]), xlab="Samples", ylab="Frequency", ylim=range(c(0,maxFreq)),  xaxt = "n", col='white')

## plot in log scale
	plot(c(log10(topPostSti[1, stimCol] + 0.0000001), log10(topPostSti[1, postStimCol] + 0.0000001)), xlab="Samples", ylab="Frequency", ylim=range(c(-7, 0)),  xaxt = "n", col='white')
	
	for (i in 1:nrow(topPostSti)) {
		points(c(log10(topPostSti[i, stimCol] + 0.0000001), log10(topPostSti[i, postStimCol] + 0.0000001)), pch=19, col=i)
		lines(c(log10(topPostSti[i, stimCol] + 0.0000001), log10(topPostSti[i, postStimCol] + 0.0000001)), lty="dotted", col=i)
		
	}

    axis(1, at=1:2, labels=c("pre-tx stim", "post-tx stim"))
    
    
	#text(1.5, 1.5, paste("npreRx: ", npreRx, "\nnpostRx: ", npostRx))
	
	dev.off()
	print("Donor reactive clones found by post-tx MLR:")
	print(paste("Total frequency in post-tx stim sample:", round(totalFreqPost, 5), sep="   "))
	print(paste("Total frequency in pre-tx stim sample:", round(totalFreqPre, 5), sep="   "))

### trend of reactive clones from pre-tx stim
	pdf(paste(prefix, "preReactive.pdf", sep = "_"))

	totalFreqPost = sum(topSti[,postStimCol])
	totalFreqPre = sum(topSti[,stimCol])
	
#	plot(c(topPostSti[1, stimCol], topPostSti[1, postStimCol]), xlab="Samples", ylab="Frequency", ylim=range(c(0,maxFreq)),  xaxt = "n", col='white')

	plot(c(log10(topSti[1, stimCol] + 0.0000001), log10(topSti[1, postStimCol] + 0.0000001)), xlab="Samples", ylab="Frequency", ylim=range(c(-7, 0)),  xaxt = "n", col='white')
	
	for (i in 1:nrow(topSti)) {
		points(c(log10(topSti[i, stimCol] + 0.0000001), log10(topSti[i, postStimCol] + 0.0000001)), pch=19, col=i)
		lines(c(log10(topSti[i, stimCol] + 0.0000001), log10(topSti[i, postStimCol] + 0.0000001)), lty="dotted", col=i)
		
	}

    axis(1, at=1:2, labels=c("pre-tx stim", "post-tx stim"))
	#text(1.5, 1.5, paste("npreRx: ", npreRx, "\nnpostRx: ", npostRx))
	dev.off()
	print("Donor reactive clones found by post-tx MLR:")
	print(paste("Total frequency in post-tx stim sample:", round(totalFreqPost, 5), sep="   "))
	print(paste("Total frequency in pre-tx stim sample:", round(totalFreqPre, 5), sep="   "))

	
### trend of reactive clones from post-tx stim in post-tx unstim samples
	pdf(paste(prefix, "unstimSamples.pdf", sep = "_"), width=10, height=10)
	par(mfrow=c(2,1))

	rename = gsub("_unstim_CD4|_unstim_CD8|_CD4|_CD8|_antidonor", "", colnames(topSti))
	rename = gsub("Subject4_|ITN4_", "S4_", rename)
	colnames(topSti) = rename

###	unstimSamples = topSti[, unstimCols]

######## Have to manually change the order !!

	newData = cbind(topSti$S4_pretx, topSti$S4_6mo, topSti$S4_12mo, topSti$S4_24mo, topSti$S4_pretx_stim, topSti$S4_12mo_posttx)
		
	totalFreq = apply(newData, 2, sum)

	newData = log10(newData + 0.0000001)
	colnames(newData) = c("pretx", "6mo", "12mo", "24mo", "pretx_stim", "12mo_posttx")


	plot(newData[1,], xlab="Samples", ylab="Frequency(log10)", col="white", xaxt="n", main="Pre-tx MLR donor reactive")

	for (i in 1:nrow(newData)) {
		points(newData[i,], pch=19, col=i)
		lines(newData[i,], lty="dotted", col=i)
	}
			
	totalFreq = apply(topSti, 2, sum)
	axis(1, at=1:ncol(newData), labels=colnames(newData))
	axis(3, at=1:ncol(newData), labels = round(totalFreq, 4))
	colnames(topPostSti) = rename

	newData = cbind(topPostSti$S4_pretx, topPostSti$S4_6mo, topPostSti$S4_12mo, topPostSti$S4_24mo, topPostSti$S4_pretx_stim, topPostSti$S4_12mo_posttx)
		
	totalFreq = apply(newData, 2, sum)

	newData = log10(newData + 0.0000001)
	colnames(newData) = c("pretx", "6mo", "12mo", "24mo", "pretx_stim", "12mo_posttx")

	
	plot(newData[1,], xlab="Samples", ylab="Frequency(log10)", col="white", xaxt="n", main="Post-tx MLR donor reactive")

	for (i in 1:nrow(newData)) {
		points(newData[i,], pch=19, col=i)
		lines(newData[i,], lty="dotted", col=i)
	}
			
	axis(1, at=1:ncol(newData), labels=colnames(newData))
	axis(3, at=1:ncol(newData), labels = round(totalFreq, 4))


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


### main function
trackClones <- function(file, freq1 = 1e-4, freq2 = 1e-5, fold = 5) {
	tcr = read.table(file, header=T)

	cd4 <- tcr[, grep("CD4|cd4", colnames(tcr))]
	cd8 <- tcr[, grep("CD8|cd8", colnames(tcr))]

	if (!is.null(cd4) & ncol(cd4) > 1) {
		cd4 <- normalize(cd4)
		compare(cd4, freq1, freq2, fold, prefix = "CD4")
	}

	if(!is.null(cd8) & ncol(cd8) > 1) { 
		cd8 <- normalize(cd8)
		compare(cd8, freq1, freq2, fold, prefix = "CD8")
	}
}
	
args<-commandArgs(TRUE)
file = args[1]
if (!is.null(file)) {
	trackClones(file)
}