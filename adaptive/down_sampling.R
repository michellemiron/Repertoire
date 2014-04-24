# downsample adaptive data

downsample <- function(tcr, target = 2000000, f = 0.5) {
	# target is the number of targeted reads per sample
	total = sum(tcr$total)
	nc = ncol(tcr)
	nsample = nc - 2
	nclone = nrow(tcr)
	
	ratio = target / total * nsample
	cat("Down-sampling to", ratio, sep="\t")
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
        "fraction" , "f", 2, "double", "[0.5] fraction of reads for T1",
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

frac = 0.5

if (!is.null(opt$fraction)) {
   frac = opt$fraction
}

if( itype == "c")  {
    tcr = read.table(file, header=T, sep=",")
} else if (itype == "t") {
  tcr = read.table(file, header=T, sep="\t")			
} else {
  tcr = read.table(file, header=T)			
}

tcr=downsample(tcr, f = frac)


write.table(tcr, "", quote=F, sep="\t", row.names = F)




