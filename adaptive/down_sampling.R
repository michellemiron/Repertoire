# downsample adaptive data

downsample <- function(tcr, target = 2000000, f = 0.5) {
	# target is the number of targeted reads per sample
#	total = sum(tcr$total)
	nc = ncol(tcr)
	nsample = nc - 2
	nclone = nrow(tcr)
	
#	ratio = target / total * nsample
#	cat("Down-sampling to", ratio, sep="\t")
	

	treads = rbinom(1, sum(tcr$total), ratio)
	readA = rmultinom(1, treads, c(f, 1-f))
	      
	# then multinomial to get frequency of each clone
	## rmultinom too slow!

	for (j in 1:2) {
	    lambda = target * tcr[,j+2] / sum(tcr[,j+2])
	    tcr[,j+2] = rpois(nclone, lambda)
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





## re-implement multinorm random number function

rmultinomial <- function(n, p){
 ncat <- length(p)
#     Check Arguments
      if (n < 0) {cat("n < 0 ","\n"); break}
      if (ncat <= 1) {cat("ncat <= 1 ","\n"); break}
      if (any(p < 0.0)) {cat("Some P(i) < 0 ","\n"); break}
      p = p / sum(p)
      if (any(p > 1.0)) {cat("Some P(i) > 1.0 ","\n"); break}
 eps <- .Machine$double.eps^0.9
      if (sum(p) > (1.0 + eps) | sum(p) < (1.0 - eps) ) {cat("Sum of P(i)
should equal 1.0 ","\n"); break}

#     Initialize variables
      ntot <- n
      sum <- 1.0
      ix <- rep(0,ncat)

#     Generate the observation
      for (icat in 1:(ncat - 1)) {
          prob <-  p[icat]/sum
          ix[icat] <-  rbinom(1,ntot,prob)
          ntot <- ntot - ix[icat]
          if (ntot <= 0) return(ix)
          sum <- sum - p[icat]
 }
      ix[ncat] <- ntot
 return (ix)
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

tcr=downsample(tcr)


write.table(tcr, "", quote=F, sep="\t", row.names = F)




