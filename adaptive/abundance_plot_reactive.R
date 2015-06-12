## do abundance plot


countUniq <- function(v) {
	b = aggregate(data.frame(count = v), list(value = v), length)
	
	return(b)

	
}

totalCounts <- function(data) {
	tc = apply(data, 2 , sum)	
	return(tc)
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

n = ncol(tcr)
mc = max(apply(tcr[,3:n],2, max))

for (i in 3:n) {
	mc[i-2] = mc[i-2] / sum(tcr[,i])
}

mf = max(mc)
fname = paste("abd.reactive", file, "pdf", sep=".")
pdf(fname, width=6, height=8)




reactive = tcr[  tcr[,3] > tcr[,5] * fold  | tcr[,4] > tcr[,6] * fold, ]
## r8 = tcr[tcr[,4] > tcr[,6] * fold, ]



# tc = totalCounts(tcr)
		
# tcr = reactive
		
flag = 0

for (i in 3:n) {
		cname = colnames(tcr)[i]
		s = reactive[,i]
		tc = sum(tcr[,i])
		# donor reactive clones

		
		counts = countUniq(s)
		if (flag == 0) {
				# plot
			plot(counts[,1]/tc, counts[,2], log="xy", col=i, xlab="clone frequency", ylab="# of clones", main=paste(file, "Donor reactive clones", sep="\n"), cex.lab=1.2, cex.axis=1.2, xlim=range(c(1e-7, 0.1)))
			flag = 1
		} else {
		points(counts[,1]/tc, counts[,2], col=i,)
		}
}
grid()
legend("topright", colnames(tcr)[3:n], col=c(3:n), pch=20)

dev.off()
