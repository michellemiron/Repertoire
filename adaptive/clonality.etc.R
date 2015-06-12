cloneCal <- function(array) {
	x = array[array >0 ] / sum(array)
#	x = sort(array, decreasing=T)
	
	l = length(x)
#	x = x[1:l]

#	x = x/ sum(x)	
#	print(x[l])
	entropy = sum(x * -1 * log2(x))
	maxentropy = -log2(1/l)
	return(signif(1 - entropy / maxentropy, 2))
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



vj = aggregate(. ~ vGeneName + jGeneName, data = tcr, sum)

#print(dim(vj))
# then calculate clonality

r1 = rep(0,4)
r2 = rep(0,4)
delta = rep(0,4)
#print(colnames(tcr)[5:8])
for (i in 5:8) {
	r1[i-4] = cloneCal(tcr[,i])
	r2[i-4] = cloneCal(vj[,i])
	delta[i-4] = r1[i-4] - r2[i-4]
}

r = cbind(colnames(tcr)[5:8], r1, r2, delta)
colnames(r) = c("samples","AA-clonality", "VJ-clonality", "delta")
write.table(t(r), quote=F)