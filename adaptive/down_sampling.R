# downsample adaptive data

downsample <- function(tcr, target = 2000000, f = 0.5) {
	# target is the number of targeted reads per sample
	total = sum(tcr$total)
	nc = ncol(tcr)
	nsample = nc - 2
	nclone = nrow(tcr)
	
	ratio = target / total * nsample
	
	# first, binomial to get the total number of reads per sample
	treads = rbinom(1, sum(tcr$total), ratio)
	readA = rmultinom(1, treads, c(f, 1-f))
	      
	# then multinomial to get frequency of each clone
	
	for (j in 1:2) {
	    tcr[,j+2] = rmultinom(1, readA[j], tcr[,j+2])
	}

	for (i in 1:nclone) {
	    tcr[i,2] = sum(tcr[i, 3:nc]) 
	}

	# normalize
	for (j in 3:nc) {
		tcr[,j] = tcr[,j] / sum(tcr[,j]) * 100
	}
	return(tcr)

}


# args<-commandArgs(TRUE)
library(getopt)

spec <- matrix(c(
        'input'     , 'i', 1, "character", "input csv or tsv file (required)",
        'type'     , 't', 2, "character", "csv or tsv (optional); default: tsv",
        'help'   , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)

opt = getopt(spec);

if (!is.null(opt$help) || is.null(opt$input)) {
    cat(paste(getopt(spec, usage=T),"\n"));
    q();
}


file = opt$input

itype = "t"
if (!is.null(opt$type)) {
   itype = opt$type
   
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

cd4 <- normalize(cd4)
cd8 <- normalize(cd8)

# first summary stats
cd4stats = repSum(cd4, fold)
cd8stats = repSum(cd8, fold)
# print(cd4stats[order(rownames(cd4stats), decreasing=T), ])
# print(cd8stats[order(rownames(cd8stats), decreasing=T), ])
# write.table(cd4stats, "", quote=F, sep="\t")
 write.table(cd8stats, "", quote=F, sep="\t")




