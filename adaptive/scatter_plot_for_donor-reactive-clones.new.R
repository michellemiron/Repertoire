## analyze TCR repertoire
# hypotheses: 
#
##1. After a certain period (6 months or later), the expanded T cell clones from donor stimulation are eliminated or reduced
##  alternatively, top clones in pre-tx stimulated condition are more likely to disappear than top clones in pre-tx unstimulated condition
##  Key caveat: sampling issue due to limited number of cells and sequencing reads. 
## solution: contorl for sampling issue by define "presence" using a threshold where the power of detection is large.

##2. The abundant clones after 6 months or later are not expanded after stimulation


compare <- function(tcr, freq1 = 1e-4, freq2 = 1e-5, fold=5, prefix = "TopClone", logscale=T) {

    if (ncol(tcr) < 2) {
        return(NA)
    }

	if (logscale == T) {
		prefix = paste(prefix, "log", sep="_")
	} 

    nsample = ncol(tcr)

	# pre-tx stim sample
    stimCol <- grep("antidonorstim|Stim|PreSt|PrSt", colnames(tcr))
    
    # pre-tx unstim sample
    preCol <- grep("pretx", colnames(tcr))[1]

	# unstim samples
	unstimCols <- grep("pretx|anti|Stim|PreSt|PrSt", colnames(tcr), invert=T)
	
     # top clones from in vitro stimulated sample
    
### pre-tx stim defined donor reactive clones; row index
    topStiRows = (tcr[,stimCol[1]] > tcr[, preCol] * fold & tcr[,stimCol[1]] >=  freq1 )
    if(length(stimCol) > 1) {
    	for (i in 2:length(stimCol)) {
    		topStiRows = (topStiRows | (tcr[,stimCol[i]] > tcr[, preCol] * fold & tcr[,stimCol[i]] >=  freq1 ))
     	}
    }
    
    topSti = tcr[topStiRows,]
    
    maxstimCol = topSti[,stimCol[1]]
    if(length(stimCol) > 1) {
		maxstimCol  = apply(topSti[,stimCol], 1, max)

	}
	
#	print(length(maxstimCol))
	rename = gsub("_cd4|_cd8", "", colnames(topSti))
	rename = gsub("ITN5|_ITN5", "", rename)

	colnames(topSti) = rename
	pretxIndex = grep("pretx", colnames(topSti))[1]
	lastIndex = grep("last", colnames(topSti))[1]
	stimIndex = grep("PreSt", colnames(topSti))[1]
	midunstimIndex = grep("pretx|anti|Stim|PreSt|PrSt|pretx|last", colnames(topSti), invert=T)
	midunstimIndex = rev(midunstimIndex)
	midunstimNames = colnames(topSti)[midunstimIndex]
	newData = cbind(maxstimCol, topSti[,pretxIndex], topSti[, midunstimIndex], topSti[,lastIndex])
	
	#	rename = gsub("Subject4_|ITN4_", "S4_", rename)
	colnames(newData) = c("pretx_stim", "pretx", midunstimNames, "last")
	
	npanel = ncol(newData) - 1
	
	pdf(paste(prefix, "frequency.pdf", sep = "_"), width=npanel * 2.9, height=4)
	par(mfrow=c(1,npanel))
	par(family="sans")
	
	
	scatterTracking(newData, freq2 = freq2, logscale = logscale)


	
	
	dev.off()
	
}

scatterTracking <- function(data, freq2 = 1e-5, logscale = T, pseudofreq = 1e-7){
	
#	data = data + 0.0000001
	maxfreq = max(data[,2:ncol(data)])

	n = nrow(data)
	for (i in 2:ncol(data)) {
		
		totalFreq = round(sum(data[,i]),6)
		
		plot(data[,1] + pseudofreq, data[,i]+ pseudofreq, log="xy", col="gray50", pch=19, xlab=colnames(data)[1], ylab="", main=paste(colnames(data)[i], "\n", "total:", totalFreq), xlim=range(c(1e-4,0.5)), ylim=range(c(pseudofreq, maxfreq)), axes = F, cex.lab=1.2)
		abline(h=freq2,col='blue', lty="dashed")

		print(freq2)			
		detectable = data[data[,i] >= freq2, ]
#		detectable = data[data[,i] >= 5e-5, ]
		n1 = nrow(detectable)
		n2 = n - n1
		
		points(detectable[,1] + pseudofreq, detectable[,i] + pseudofreq, col="salmon", pch=19)
#		text(1e-1, 1e-4, labels = paste("N =", n1), col='red')
#		text(1e-1, 1e-6, labels = paste("N =", n2), col='gray30')

		axis(1,at=c(1,10,100, 1000)*1e-4, cex.axis=1.3, cex.names= 1.2)
		axis(2, at = c(1,10, 100, 1000, 10000)* pseudofreq, labels=c("0", c(10, 100, 1000, 10000) * pseudofreq), cex.axis=1.3, cex.names= 1.2)
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
	
	plot(newData[1,], xlab="Samples", ylab=ylabel, col="white", xaxt="n", main=title, ylim = yrange)

	for (i in 1:nrow(newData)) {
		points(newData[i,], pch=19, col=i)
		lines(newData[i,], lty="dotted", col=i)
	}
			
	axis(1, at=1:ncol(newData), labels=colnames(newData))
	axis(3, at=1:ncol(newData), labels = round(totalFreq, 4), tick = F, outer =F, tcl=0)
	
	
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
	sample = unlist(strsplit(file, "_"))[1]
	trackClones(tcr,freq1 = freq1, freq2 = freq2, fold, prefix = sample,  logscale = logscale)
}
	
trackClones <- function(tcr, freq1 = 1e-4, freq2 = 1e-5, fold = 5, prefix = "Sample", logscale=T) {
	cd4 <- tcr[, grep("CD4|cd4", colnames(tcr))]
	cd8 <- tcr[, grep("CD8|cd8", colnames(tcr))]

	if (!is.null(cd4) & ncol(cd4) > 1) {
		cd4 <- normalize(cd4)
		compare(cd4, freq1 = freq1, freq2 = freq2, fold = fold, prefix = paste(prefix, "CD4", sep="_"), logscale)
	}

	if(!is.null(cd8) & ncol(cd8) > 1) { 
		cd8 <- normalize(cd8)
		compare(cd8, freq1 = freq1, freq2 = freq2, fold = fold, prefix = paste(prefix, "CD8", sep="_"), logscale)
	}

	
}
	
args<-commandArgs(TRUE)
file = args[1]
freq2 = args[2]

if (is.null(freq2)) {
	freq2 = 1e-5
} else {
	freq2 = as.numeric(freq2)
}		
trackClonesCL(file, freq2 = freq2, logscale = T)
