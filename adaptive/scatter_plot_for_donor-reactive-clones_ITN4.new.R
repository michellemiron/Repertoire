## analyze TCR repertoire
# hypotheses: 
#
##1. After a certain period (6 months or later), the expanded T cell clones from donor stimulation are eliminated or reduced
##  alternatively, top clones in pre-tx stimulated condition are more likely to disappear than top clones in pre-tx unstimulated condition
##  Key caveat: sampling issue due to limited number of cells and sequencing reads. 
## solution: contorl for sampling issue by define "presence" using a threshold where the power of detection is large.

##2. The abundant clones after 6 months or later are not expanded after stimulation


compare <- function(tcr, freq1, freq2, fold=5, prefix = "TopClone", logscale=T) {
    if (ncol(tcr) < 2) {
        return(NA)
    }

	if (logscale == T) {
		prefix = paste(prefix, "log", sep="_")
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



	rename = gsub("_unstim_CD4|_unstim_CD8|_CD4|_CD8|_antidonor", "", colnames(topSti))
	rename = gsub("Subject4_|ITN4_", "S4_", rename)
	
### trend of reactive clones from post-tx stim in post-tx unstim samples
	pdf(paste(prefix, "frequency.pdf", sep = "_"), width=11.6 , height=4)
	par(mfrow=c(1,4))
	par(family="sans")
	
######## Have to manually change the order !!
	
## panel 1, freq tracking for donor reactive clones from pre-tx MLR

	colnames(topSti) = rename
	newData = cbind(topSti$S4_pretx_stim, topSti$S4_pretx, topSti$S4_6mo, topSti$S4_12mo, topSti$S4_24mo)
	colnames(newData) = c("pretx_stim", "pretx", "6mo", "12mo", "24mo")
	scatterTracking(newData, freq2, logscale = logscale, title = "Pre-tx MLR donor reactive clones\n")


## panel 2, freq tracking for donor reactive clones from post-tx MLR
#	colnames(topPostSti) = rename
#	newData = cbind( topPostSti$S4_12mo_posttx, topPostSti$S4_pretx, topPostSti$S4_6mo, topPostSti$S4_12mo, topPostSti$S4_24mo, topPostSti$S4_pretx_stim)
#	colnames(newData) = c("12mo_posttx", "pretx", "6mo", "12mo", "24mo", "pretx_stim" )
#	scatterTracking(newData, logscale = logscale, title = "Post-tx MLR donor reactive clones\n")
	
	
	dev.off()
	
}

scatterTracking <- function(data, freq2 = 1e-5, title = "", logscale = T, pseudofreq = 1e-7){
	
#	data = data + 0.0000001
	maxfreq = max(data[,2:ncol(data)])
	n = nrow(data)
	for (i in 2:ncol(data)) {
		detectable = data[data[,i] >= freq2, ]
		totalFreq = round(sum(data[,i]),4)
		n1 = nrow(data[data[,i] >= freq2, ])
		n2 = n - n1
		plot(data[,1] + pseudofreq, data[,i]+ pseudofreq, log="xy", col="gray50", pch=19, xlab=colnames(data)[1], ylab="", main=paste(colnames(data)[i], "\n", "total:", totalFreq), xlim=range(c(1e-4,0.5)), ylim=range(c(pseudofreq, maxfreq)), axes = F, cex.lab=1.2)
		abline(h=freq2,col='blue', lty="dashed")
		points(detectable[,1] + pseudofreq, detectable[,i] + pseudofreq, col="salmon", pch=19)
#		text(1e-1, 1e-4, labels = paste("N =", n1), col='red')
#		text(1e-1, 1e-6, labels = paste("N =", n2), col='gray30')

		axis(1,at=c(1,10,100, 1000)*1e-4, cex.axis=1.3,)
		axis(2, at = c(1,10, 100, 1000, 10000)* pseudofreq, labels=c("0", c(10, 100, 1000, 10000) * pseudofreq), cex.axis=1.3)
	}
}


plotTracking <- function(newData, title = "", logscale = T) {
	zeroadd = 0.00000001 
	totalFreq = apply(newData, 2, sum)
	ylabel = "Frequency"
	yrange = range(c(0, max(newData)))

	if (logscale == T) {
		newData = log10(newData + zeroadd)
		ylabel = "Frequency(log10)"
		yrange = range(c(log10(zeroadd), 0))
	} 
	
	plot(newData[1,], xlab="Samples", ylab=ylabel, col="white", xaxt="n", main=title, ylim = yrange, cex.lab=1.2)

	for (i in 1:nrow(newData)) {
		points(newData[i,], pch=19, col=i)
		lines(newData[i,], lty="dotted", col=i)
	}
			
	axis(1, at=1:ncol(newData), labels=colnames(newData), cex.axis=1.5)
	axis(3, at=1:ncol(newData), labels = round(totalFreq, 4), tick = F, outer =F, tcl=0, cex.axis=1.5)
	
	
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
trackClonesCL <- function(file, freq1 = 1e-4, freq2 = 1e-5, fold = 5, logscale = T) {
	tcr = read.table(file, header=T, check.names=FALSE)
	trackClones(tcr,freq1, freq2, fold, logscale = logscale)
}
	
trackClones <- function(tcr, freq1 = 1e-4, freq2 = 1e-5, fold = 5, logscale=T) {
	cd4 <- tcr[, grep("CD4|cd4", colnames(tcr))]
	cd8 <- tcr[, grep("CD8|cd8", colnames(tcr))]

	if (!is.null(cd4) & ncol(cd4) > 1) {
		cd4 <- normalize(cd4)
		compare(cd4, freq1, freq2, fold, prefix = "CD4", logscale)
	}

	if(!is.null(cd8) & ncol(cd8) > 1) { 
		cd8 <- normalize(cd8)
		compare(cd8, freq1, freq2, fold, prefix = "CD8", logscale)
	}

	
}
	
args<-commandArgs(TRUE)
file = args[1]
logscale = args[2]

if (!is.null(file)) {
	if ( !is.null(logscale)) {
		logflag = T
	} else {
		logflag = F
	}
		
	trackClonesCL(file, logscale = logflag)
}