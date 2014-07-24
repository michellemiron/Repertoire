## analyze TCR repertoire
# hypotheses: 
#
##1. After a certain period (6 months or later), the expanded T cell clones from donor stimulation are eliminated or reduced
##  alternatively, top clones in pre-tx stimulated condition are more likely to disappear than top clones in pre-tx unstimulated condition
##  Key caveat: sampling issue due to limited number of cells and sequencing reads. 
## solution: contorl for sampling issue by define "presence" using a threshold where the power of detection is large.

##2. The abundant clones after 6 months or later are not expanded after stimulation

r50Cal <- function(array) {
	x = array[array >0 ] / sum(array)
	s = sort(x, decreasing=T)
	l = length(x)
	total = 0
	for (i in 1:l) {
		total = total + s[i]
		if (total >= 0.5) {
			return(signif(i/l,2))
			
		}
	}
	return(1.0)
}

r20Cal <- function(array) {
	x = array[array >0 ] / sum(array)
	s = sort(x, decreasing=T)
	l = length(x)
	total = 0
	for (i in 1:l) {
		total = total + s[i]
		if (total >= 0.2) {
			return(signif(i/l,2))
			
		}
	}
	return(1.0)
}


# clonality
cloneCal <- function(array) {
	x = array[array >0 ] / sum(array)
	l = length(x)
	total = 0
	entropy = sum(x * -1 * log2(x))
	maxentropy = -log2(1/l)
	return(signif(1 - entropy / maxentropy, 2))
}

# summary statistics

repSum <- function(tcr) {
       if (ncol(tcr) < 1) {
       return(NA)
}   else {   
	cols = colnames(tcr)
	
	r  = matrix(nrow = length(cols), ncol = 6)
	rownames(r) = cols
	colnames(r) = c("N_clones", "Min_freq", "Max_freq", "R50", "R20", "Clonality")
	
	for  (i in 1:length(cols)) {
		coln = cols[i]
		nclones = length(tcr[tcr[,i]>0, i])
		minfreq = min(tcr[tcr[,i]>0, i])/sum(tcr[,i])
		maxfreq = signif(max(tcr[,i])/sum(tcr[,i]), 3)

		r50 = r50Cal(tcr[,i])
		r20 = r20Cal(tcr[,i])
		clonality = cloneCal(tcr[,i])
		r[i,] = c(nclones, minfreq, maxfreq, r50, r20, clonality)
	}
	
	return(r)
}
}

shannon.entropy <- function(p)
{
	if (min(p) < 0 || sum(p) <= 0)
		return(NA)
	p.norm <- p[p>0]/sum(p)
	-sum(log2(p.norm)*p.norm)
}


jensen_shannon <- function(p, q){
	## JSD = H(0.5 *(p +q)) - 0.5H(p) - 0.5H(q)
 	# H(X) = \sum x_i * log2(x_i)
#	p = p[p >0 & q >0]
#	q = q[p>0 & q>0]
	p = p / sum(p)
	q = q / sum(q)
	Hj = shannon.entropy(0.5 *(p+q)) 
	Hp = shannon.entropy(p) 
	Hq = shannon.entropy(q)
	
	jsd = Hj - 0.5*(Hp+Hq)
#	cat(Hj, Hp, Hq, jsd, "\n")
	return(jsd)
}


compare <- function(tcr, freq1, freq2, fold=10) {
#	print(dim(tcr))
	if (ncol(tcr) < 2) {
		return(NA)
	}
#	x = tcr[which(rowMeans(tcr) > 0.001),]
	
	## need to renormalize!! 
	x = normalize(tcr)
	# x = tcr
#	print(dim(x))
	nsample = ncol(x)
#	jsd <- matrix(ncol=nsample, nrow=nsample)
#	colnames(jsd) = colnames(tcr)
#	rownames(jsd) = colnames(tcr)
	##pair-wise jensen-shannon divergence
#	for (i in 1:(nsample-1)){
	#	cat(colnames(tcr)[i], "\n")	
#		jsd[i,i] = 0
#		for (j in (i+1):nsample){
#			jsd[j,j] = 0
#			jsd[i,j] = round(jensen_shannon(x[,i]/100, x[,j]/100), 3)
#			jsd[j,i] = jsd[i,j]
#		}
#	}
#	cat("J-S D on CDR3", "\n")
#	print(jsd)
		
#	cat(colnames(tcr), "\n")

#        stimCol <- grep("stim|PreSt|PrSt|Stim", colnames(tcr))[1]
        stimCol <- grep("stim|PreSt|PrSt|Stim", colnames(tcr))
        
        preCol <- grep("tx|pre", colnames(tcr))[1]
        postCols <- grep("stim|pre|tx|PreSt|PrSt|Stim", colnames(tcr), invert=T)


        

# top clones from in vitro stimulated sample
        topStiRows = (tcr[,stimCol[1]] > tcr[, preCol] * fold & tcr[,stimCol[1]] >= 100 * freq1 )
#        topStiRows = (tcr[,stimCol[2]] > tcr[, preCol] * fold & tcr[,stimCol[2]] >= 100 * freq1 )	

        if (length(stimCol) > 1) {
            for (i in 2:length(stimCol)) {
                topStiRows = (topStiRows | (tcr[,stimCol[i]] > tcr[, preCol] * fold & tcr[,stimCol[i]] >= 100 * freq1 ))
            }
        }
        
        # topSti = tcr[tcr[,stimCol] > tcr[, preCol] * fold & tcr[,stimCol] >= 100 * freq1, ]
        topSti = tcr[topStiRows,]
                                        #	cat(dim(topSti), "\n") 
	
#	cat(dim(top3pSti), "\n")
	

## test 1: look at top expanded clones in stimulated condition, test if they are reduced in post
#	cat("Test if top expanded clones are reduced even than pre-tx unsti\n")

	# preFreq = topSti[topSti[,preCol] >= freq2 * 100, preCol]


	pretx = c(length(which(topSti[, preCol ] >= freq2 * 100)), length(which(topSti[, preCol] < freq2 * 100)))

        
	cat(colnames(tcr)[preCol],  pretx, "\n", sep="\t")

#	for (i in postCols) {
#		y = c(length(which(topSti[, i] >= freq2 * 100)), length(which(topSti[, i] < freq2 * 100 )))
#		test = fisher.test(cbind(y, pretx))
#		cat(colnames(tcr)[i], y,signif(test$p.value,2),  round(test$estimate,2), "\n", sep="\t")
#	}

#	cat("------------", "\n")
## test 2: 3rd-party sti
	tpStimCol <- grep("3rd_party", colnames(tcr))

#	cat(tpStimCol, "\n")
	if (length(tpStimCol) > 0 ) {
            tpStimCol <- tpStimCol[1]

            top3pSti = tcr[tcr[,tpStimCol] > tcr[,preCol] * fold & tcr[, tpStimCol] >= 100 * freq1, ]

           # preFreq3p = top3pSti[top3pSti[,preCol] >= freq2 * 100, preCol]		
            pretx3p = c(length(which(top3pSti[, preCol ] >= freq2 * 100)), length(which(top3pSti[, preCol] < freq2 * 100)))

#		cat("Third-party stimulation\n")
#		cat(colnames(tcr)[preCol], " tp:", pretx3p, "\n")

            for (i in postCols) {
                y3p = c(length(which(top3pSti[, i] >= freq2 * 100)), length(which(top3pSti[, i] < freq2 * 100 )))
                test3p = fisher.test(cbind(y3p, pretx3p))
                cat(colnames(tcr)[i], y3p,  signif(test3p$p.value,2),  round(test3p$estimate,2), "\n", sep="\t")
                
            }
        }
        cat("------------", "\n")
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

#if (opt$vj != 0) {
#	vj = aggregate(. ~ vGeneName + jGeneName, data = tcr, sum)

#	vjDiverg(vj[,grep("CD4", colnames(tcr))])
#	vjDiverg(vj[,grep("CD8", colnames(tcr))])
#}

cd4 <- tcr[, grep("CD4|cd4", colnames(tcr))]
cd8 <- tcr[, grep("CD8|cd8", colnames(tcr))]

cd4 <- normalize(cd4)
cd8 <- normalize(cd8)

# first summary stats
#cd4stats = repSum(cd4)
cd8stats = repSum(cd8)

#print(cd4stats[order(rownames(cd4stats), decreasing=T), ])
print(cd8stats[order(rownames(cd8stats), decreasing=T), ])

if (!is.null(cd4) & ncol(cd4) > 1) {
compare(cd4, freq1, freq2, fold)
}

if(!is.null(cd8) & ncol(cd8) > 1) { 
#compare(cd8, freq1, freq2, fold)
}

## compute distance



# quit(save = "no")

