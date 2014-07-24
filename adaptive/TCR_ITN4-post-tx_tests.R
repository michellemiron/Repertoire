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

    postCols <- grep("pretx", colnames(tcr), invert=T)
    # postCols <- grep("tim|tx", colnames(tcr), invert=T)
    
                                        # top clones from in vitro stimulated sample
    # topSti = tcr[tcr[,stimCol] > tcr[, preCol] * fold & tcr[,stimCol] >= 100 * freq1, ]

    topStiRows = (tcr[,stimCol[1]] > tcr[, preCol] * fold & tcr[,stimCol[1]] >= 100 * freq1 )
    if (length(stimCol) > 1) {
        for (i in 2:length(stimCol)) {
            topStiRows = (topStiRows | (tcr[,stimCol[i]] > tcr[, preCol] * fold & tcr[,stimCol[i]] >= 100 * freq1 ))
        }
    }
    
                                        # topSti = tcr[tcr[,stimCol] > tcr[, preCol] * fold & tcr[,stimCol] >= 100 * freq1, ]
    topSti = tcr[topStiRows,]
    allPre = tcr[tcr[,preCol] >= freq2 * 100, ]
    
	

## test 1: look at top expanded clones in stimulated condition, test if they are reduced in post
#	cat("Test if top expanded clones are reduced even than pre-tx unsti\n")

	# preFreq = topSti[topSti[,preCol] >= freq2 * 100, preCol]


	pretx = c(length(which(topSti[, preCol ] >= freq2 * 100)), length(which(topSti[, preCol] < freq2 * 100)))

 	cat(colnames(tcr)[preCol], ":", pretx, "\n")
#	cat(mean(preFreq), "\n")
    for (i in postCols) {
        
#        x = c(length(which(allPre[, i] >= freq2 * 100)), length(which(allPre[, i] < freq2 * 100 )))
        y = c(length(which(topSti[, i] >= freq2 * 100)), length(which(topSti[, i] < freq2 * 100 )))
        test = fisher.test(cbind(y, pretx))
        cat(colnames(tcr)[i], ":", pretx, y, signif(test$p.value,2),  round(test$estimate,2), round(test$conf.int[1:2], 2),  "\n")
        
	
    }
    
}


normalize <- function(tcr){
    if (ncol(tcr) < 1) {
        return(NA)
    }
    for (f in 1:ncol(tcr)) {
		tcr[,f] = tcr[,f] / sum(tcr[,f]) * 100

	}
	return(tcr)	
}


vjDiverg <- function(vj) {
	if (ncol(vj) < 1) {
		return(NA)
	}
	
#	x = tcr[which(rowMeans() > 0.001),]
	
	## need to renormalize!! 
	x = normalize(vj)
	# x = tcr
	print(dim(x))
	x.sum <- repSum(x)
	print(x.sum[order(rownames(x.sum), decreasing=T),])
	x = x[which(rowMeans(x)> 0.001), ]
	
	nsample = ncol(x)
	jsd <- matrix(ncol=nsample, nrow=nsample)
	colnames(jsd) = colnames(x)
	rownames(jsd) = colnames(x)
	##pair-wise jensen-shannon divergence
	for (i in 1:(nsample-1)){
	#	cat(colnames(x)[i], "\n")	
		jsd[i,i] = 0
		for (j in (i+1):nsample){
			jsd[j,j] = 0
			jsd[i,j] = round(jensen_shannon(x[,i]/100, x[,j]/100), 3)
			jsd[j,i] = jsd[i,j]
		}
	}

	print(jsd)

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

