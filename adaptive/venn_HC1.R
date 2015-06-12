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
    
    T1STIMCol <- grep("T1_STIM", colnames(tcr))
    T1UNCol <- grep("T1_UN", colnames(tcr))
	
	T2STIMCol <- grep("T2_STIM", colnames(tcr))
	T2UNCol = T1UNCol
    

    T1topStiRows = (tcr[, T1STIMCol] > tcr[, T1UNCol] * fold & tcr[, T1STIMCol] >= 100 * freq1 )
    T2topStiRows = (tcr[, T2STIMCol] > tcr[, T1UNCol] * fold & tcr[, T2STIMCol] >= 100 * freq1 )
    
    
	x = nrow(tcr[T1topStiRows, ])
	y = nrow(tcr[T2topStiRows, ])
	z = nrow(tcr[T1topStiRows & T2topStiRows, ])
	
	freq1 = sum(tcr[T1topStiRows, T1STIMCol])
	freq2 = sum(tcr[T2topStiRows, T2STIMCol])
	freqO1 = sum(tcr[T1topStiRows  & T2topStiRows, T1STIMCol])
	freqO2 = sum(tcr[T1topStiRows  & T2topStiRows, T2STIMCol])
	
	freq1p = sum(tcr[T1topStiRows, T1UNCol])
	freq2p = sum(tcr[T2topStiRows, T1UNCol])
	
	
	cat(x,y,z, "\n", sep="\t")
	cat(freq1, freq2, freqO1, freqO2, freq1p, freq2p, "\n", sep="\t")    

	   
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
fold = 5  # fold change default threshold

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


compare(cd4, freq1, freq2, fold)

compare(cd8, freq1, freq2, fold)

