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
	
	## need to renormalize!! 
	x = normalize(tcr)

	nsample = ncol(x)

	stimCol <- grep("stim_CD|stim2|stim3", colnames(tcr))
        
        preCol <- grep("pre", colnames(tcr))[1]
        postCols <- grep("pre|stim|anti", colnames(tcr), invert=T)
       

# top clones from in vitro stimulated sample
        topStiRows = (tcr[,stimCol[1]] > tcr[, preCol] * fold & tcr[,stimCol[1]] >= 100 * freq1 )

        if (length(stimCol) > 1) {
            for (i in 2:length(stimCol)) {
                topStiRows = (topStiRows | (tcr[,stimCol[i]] > tcr[, preCol] * fold & tcr[,stimCol[i]] >= 100 * freq1 ))
            }
        }
        
        topSti = tcr[topStiRows,]

	tpStimCol <- grep("3rd_party", colnames(tcr))
        top3pSti = tcr[tcr[,tpStimCol] > tcr[,preCol] * fold & tcr[, tpStimCol] >= 100 * freq1, ]

        pretx3p = c(length(which(top3pSti[, preCol ] >= freq2 * 100)), length(which(top3pSti[, preCol] < freq2 * 100)))

	pretx = c(length(which(topSti[, preCol ] >= freq2 * 100)), length(which(topSti[, preCol] < freq2 * 100)))
	cat(colnames(tcr)[preCol],  pretx, "", "",pretx3p, "\n", sep="\t")

	for (i in postCols) {
		y = c(length(which(topSti[, i] >= freq2 * 100)), length(which(topSti[, i] < freq2 * 100 )))
		test = fisher.test(cbind(y, pretx))
                y3p = c(length(which(top3pSti[, i] >= freq2 * 100)), length(which(top3pSti[, i] < freq2 * 100 )))
                test3p = fisher.test(cbind(y3p, pretx3p))
#		ctest = fisher.test(cbind(y,z))	    
		cat(colnames(tcr)[i], y,signif(test$p.value,2),  round(test$estimate,2), y3p,  signif(test3p$p.value,2),  round(test3p$estimate,2),  "\n", sep="\t")
	}


        cat("------------", "\n")
}


normalize <- function(tcr){
	if (ncol(tcr) < 1) {
		return(NA)
	}
#	cat(colnames(tcr), "\n")
#	cat(c(tcr[1:3, ]), "\n")

	for (f in 1:ncol(tcr)) {
		tcr[,f] = tcr[,f] / sum(tcr[,f], na.rm = T) * 100
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


cd4 <- tcr[, grep("CD4|cd4", colnames(tcr))]
cd8 <- tcr[, grep("CD8|cd8", colnames(tcr))]

#cd4 <- normalize(cd4)
#cd8 <- normalize(cd8)

compare(cd4, freq1, freq2, fold)

compare(cd8, freq1, freq2, fold)
